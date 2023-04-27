import pybamm
import numpy as np
import matplotlib.pylab as plt
import pints
import time


class SingleReactionSolution(pints.ForwardModel):
    def __init__(self, n=100, max_x=10, param=None):
        # Set fixed parameters here
        if param is None:
            param = pybamm.ParameterValues({
                "Faraday Constant [C mol-1]": 96485.3328959,
                "Gas constant [J K-1 mol-1]": 8.314459848,
                "Far-field concentration of A [mol cm-3]": 2e-6,
                "Diffusion Constant [cm2 s-1]": 1e-5,
                "Electrode Area [cm2]": 1,
                "Temperature [K]": 298.2,
                "Voltage frequency [rad s-1]": 9.0152,
                "Voltage start [V]": 0.5,
                "Voltage reverse [V]": -0.5,
                "Voltage amplitude [V]": 0.0,
                "Scan Rate [V s-1]": 0.05,
                "Electrode Coverage [mol cm-2]": 1e-12,
            })

        # Create dimensional fixed parameters
        c_inf = pybamm.Parameter("Far-field concentration of A [mol cm-3]")
        D = pybamm.Parameter("Diffusion Constant [cm2 s-1]")
        F = pybamm.Parameter("Faraday Constant [C mol-1]")
        R = pybamm.Parameter("Gas constant [J K-1 mol-1]")
        S = pybamm.Parameter("Electrode Area [cm2]")
        T = pybamm.Parameter("Temperature [K]")

        E_start_d = pybamm.Parameter("Voltage start [V]")
        E_reverse_d = pybamm.Parameter("Voltage reverse [V]")
        #deltaE_d = pybamm.Parameter("Voltage amplitude [V]")
        v = pybamm.Parameter("Scan Rate [V s-1]")

        # Create dimensional input parameters
        E0 = pybamm.InputParameter("Reversible Potential [non-dim]")
        k0 = pybamm.InputParameter("Reaction Rate [non-dim]")
        alpha = pybamm.InputParameter("Symmetry factor [non-dim]")
        Cdl = pybamm.InputParameter("Capacitance [non-dim]")
        Ru = pybamm.InputParameter("Uncompensated Resistance [non-dim]")
        #omega_d = pybamm.InputParameter("Voltage frequency [rad s-1]")


        E0_d = pybamm.InputParameter("Reversible Potential [V]")
        k0_d = pybamm.InputParameter("Reaction Rate [s-1]")
        alpha = pybamm.InputParameter("Symmetry factor [non-dim]")
        Cdl_d = pybamm.InputParameter("Capacitance [F]")
        Ru_d = pybamm.InputParameter("Uncompensated Resistance [Ohm]")

        # Create scaling factors for non-dimensionalisation
        E_0 = R * T / F
        T_0 = E_0 / v
        L_0 = pybamm.sqrt(D * T_0)
        I_0 = D * F * S * c_inf / L_0

        # Non-dimensionalise parameters
        E0 = E0_d / E_0
        k0 = (k0_d * pybamm.sqrt(S))/D
        Cdl = Cdl_d * S * E_0 / (I_0 * T_0)
        Ru = Ru_d * I_0 / E_0
        #omega = 2 * np.pi * omega_d * T_0

        E_start = E_start_d / E_0
        E_reverse = E_reverse_d / E_0
        t_reverse = E_start - E_reverse
        #deltaE = deltaE_d / E_0

        # Input voltage protocol
        Edc_forward = -pybamm.t
        Edc_backwards = pybamm.t - 2*t_reverse
        Eapp = E_start + \
            (pybamm.t <= t_reverse) * Edc_forward + \
            (pybamm.t > t_reverse) * Edc_backwards 
            #deltaE * pybamm.sin(omega * pybamm.t)

        # create PyBaMM model object
        model = pybamm.BaseModel()

        # Create state variables for model
        theta = pybamm.Variable("ratio_A", domain="solution")
        i = pybamm.Variable("Current")

        # Effective potential
        Eeff = Eapp - i * Ru

        # Faradaic current
        i_f = pybamm.BoundaryGradient(theta, "left")

        # ODE equations
        model.rhs = {
            theta: pybamm.div(pybamm.grad(theta)),
            i: (i_f + Cdl * Eapp.diff(pybamm.t) - i), # 1/(Cdl * Ru) *
        }

        # algebraic equations (none)
        model.algebraic = {
        }

        # Butler-volmer boundary condition at electrode
        theta_at_electrode = pybamm.BoundaryValue(theta, "left")
        butler_volmer = k0 * pybamm.exp(-alpha * (Eeff - E0)) * (theta_at_electrode 
                        - ((1 - theta_at_electrode)*pybamm.exp(Eeff - E0)))

        # Boundary and initial conditions
        model.boundary_conditions = {
            theta: {
                "right": (pybamm.Scalar(1), "Dirichlet"),
                "left": (butler_volmer, "Neumann"),
            }
        }

        model.initial_conditions = {
            theta: pybamm.Scalar(1),
            i: Cdl * (1.0),
        }

        # set spatial variables and solution domain geometry
        x = pybamm.SpatialVariable('x', domain="solution")
        default_geometry = pybamm.Geometry({
            "solution": {
                x: {"min": pybamm.Scalar(0), "max": pybamm.Scalar(max_x)}
            }
        })

        default_var_pts = {
            x: n
        }

        # Using Finite Volume discretisation on an expanding 1D grid for solution
        default_submesh_types = {
            "solution": pybamm.MeshGenerator(pybamm.Exponential1DSubMesh, {'side': 'left'})
        }
        default_spatial_methods = {
            "solution": pybamm.FiniteVolume()
        }

        # model variables
        model.variables = {
            "Current [non-dim]": i,
            "Applied Voltage [non-dim]": Eapp,
        }

        #--------------------------------

        # Set model parameters
        param.process_model(model)
        geometry = default_geometry
        param.process_geometry(geometry)

        # Create mesh and discretise model
        mesh = pybamm.Mesh(geometry, default_submesh_types, default_var_pts)
        disc = pybamm.Discretisation(mesh, default_spatial_methods)
        disc.process_model(model)

        # Create solver
        #solver = pybamm.CasadiSolver(mode="fast",
        #                             rtol=1e-9,
        #                             atol=1e-9,
        #                             extra_options_setup={'print_stats': False})
        #model.convert_to_format = 'jax'
        #solver = pybamm.JaxSolver(method='BDF')
        #model.convert_to_format = 'python'
        #solver = pybamm.ScipySolver(method='BDF')
        
        model.convert_to_format = 'python'
        solver = pybamm.ScipySolver(method='Radau', rtol=1e-6, atol=1e-6)

        # Store discretised model and solver
        self._model = model
        self._solver = solver
        self._fast_solver = None
        self._omega_d = param["Voltage frequency [rad s-1]"]

        self._I_0 = param.process_symbol(I_0).evaluate()
        self._T_0 = param.process_symbol(T_0).evaluate()
        self._E_0 = param.process_symbol(E_0).evaluate()
        self._L_0 = param.process_symbol(L_0).evaluate()
        self._S = param.process_symbol(S).evaluate()
        self._D = param.process_symbol(D).evaluate()
        self._default_var_points = default_var_pts

    def non_dim(self, x):
        if len(x) > 6:
            return np.array([
                x[0] * self._L_0 / self._D,
                x[1] / self._E_0,
                x[2],
                x[3] * self._I_0 / self._E_0,
                x[4] * self._S * self._E_0 / (self._I_0 * self._T_0),
                x[5] * 2 * np.pi * self._T_0,
                x[6]
            ])
        else:
            return np.array([
                x[0] * self._L_0 / self._D,
                x[1]/self._E_0,
                x[2],
                x[3] * self._I_0 / self._E_0,
                x[4] * self._S * self._E_0 / (self._I_0 * self._T_0),
                x[5] * 2 * np.pi * self._T_0
            ])

    def simulate(self, parameters, times):
        input_parameters = {
            "Reaction Rate [s-1]": parameters[0],
            "Reversible Potential [V]": parameters[1],
            "Symmetry factor [non-dim]": parameters[2],
            "Uncompensated Resistance [Ohm]": parameters[3],
            "Capacitance [F]": parameters[4],
            "Voltage frequency [rad s-1]": parameters[5],
        }
        # if self.fast_solver is None:
        #    solution = self._solver.solve(self._model, times, inputs=input_parameters)

        index = list(self._default_var_points.values())[0]

        try:
            solution = self._solver.solve(self._model, times, inputs=input_parameters)
            return (
                solution["Current [non-dim]"](times),
                solution["Applied Voltage [non-dim]"](times),
                #parameters
            )
        except pybamm.SolverError:
            print('solver errored for params',parameters)
            return np.zeros_like(times)

    def n_parameters(self):
        return 6


if __name__ == '__main__':
    

    #FOR DIGIELCHM comparison
# =============================================================================
#     files = [['digielchcomp/current k0 1e-3.txt',24234345325654264564.0], ['digielchcomp/current k0 2e-3.txt', 2],
#               ['digielchcomp/current k0 5e-3.txt', 5], ['digielchcomp/current k0_0.1.txt', 100],
#               ['digielchcomp/current k0_0.5.txt', 500], ['digielchcomp/current k0_1.txt', 1000],
#               ['digielchcomp/current k0_1e-2.txt', 10], ['digielchcomp/current k0_10.txt', 10000],
#               ['digielchcomp/current k0_100.txt', 100000], ['digielchcomp/current k0_1000.txt', 1000000]]
#     for i in files:
# =============================================================================
        
# =============================================================================
#         print(i)
#         
#         volt = []
#         curr = []
#         row = []
#           
#         f = open(i[0],'r')
#         for row in f:
#             row = row.split("\t")
#             volt.append(float(row[0]))
#             curr.append(float(row[1]))
# =============================================================================
            
            
        # pybamm.set_logging_level('INFO')
        model = SingleReactionSolution()
        
        React_Rate = 1
        Revers_Potent = 0
        Symm_fact = 0.5
        Uncomp_resist = 1.0
        Capacit = 1e-8
        Volt_freq = 0
        x = np.array([React_Rate, Revers_Potent, Symm_fact, Uncomp_resist, Capacit, Volt_freq])
        
        Tmax = abs(0.5 + 0.5)/0.05 * 2
        Tdim = np.linspace(0, Tmax, 2**8)
        TnonDim = (96485.3328959 * 0.05 / (8.314459848*298)) * Tdim

        t0 = time.perf_counter()
        current, voltage = model.simulate(x, TnonDim)
        t1 = time.perf_counter()
                
        ##redimensionalizing here for now. Messy to do in main, move later
        I_d = current*model._I_0*1e2
        E_d = voltage * model._E_0
        #title = i[0]
          
        
        plt.cla()
        #plt.plot(volt, curr, color = 'g', label = 'Digielch')
        plt.plot(E_d, I_d, color = 'b', label = 'Pybamm')
        #plt.title(title)
        plt.legend()
        plt.xlabel("Eapp [V]")
        plt.ylabel("current [A]")
        plt.show()
        #plt.savefig(i[0]+".png")

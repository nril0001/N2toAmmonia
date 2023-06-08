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
                "Far-field concentration of A [mol cm-3]": 1e-6,
                "Diffusion Constant [cm2 s-1]": 600,
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
        v = pybamm.Parameter("Scan Rate [V s-1]")

        # Create dimensional input parameters
        E0 = pybamm.InputParameter("Reversible Potential [non-dim]")
        k0 = pybamm.InputParameter("Reaction Rate [non-dim]")
        alpha = pybamm.InputParameter("Symmetry factor [non-dim]")
        Cdl = pybamm.InputParameter("Capacitance [non-dim]")
        Ru = pybamm.InputParameter("Uncompensated Resistance [non-dim]")


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

        E_start = E_start_d / E_0
        E_reverse = E_reverse_d / E_0
        t_reverse = E_start - E_reverse
        
        #creating time scale and non-dimensionalizing
        Tmax_nd = ((abs(E_start_d - E_reverse_d) / v * 2)/ T_0)
        
        #number of time steps
        m = 20000
        
        #length of time step, nondimensional
        deltaT_nd = Tmax_nd / m

        # Input voltage protocol
        Edc_forward = -pybamm.t
        Edc_backwards = pybamm.t - 2*t_reverse
        Eapp = E_start + \
            (pybamm.t <= t_reverse) * Edc_forward + \
            (pybamm.t > t_reverse) * Edc_backwards 

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
            theta: pybamm.Scalar(0),
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
            "theta": theta_at_electrode,
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

        self._I_0 = param.process_symbol(I_0).evaluate()
        self._T_0 = param.process_symbol(T_0).evaluate()
        self._E_0 = param.process_symbol(E_0).evaluate()
        self._L_0 = param.process_symbol(L_0).evaluate()
        self._S = param.process_symbol(S).evaluate()
        self._D = param.process_symbol(D).evaluate()
        self._default_var_points = default_var_pts
        self._deltaT_nd = param.process_symbol(deltaT_nd).evaluate()

        # store time scale related things
        self._Tmax_nd = param.process_symbol(Tmax_nd).evaluate()
        self._m = m

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
                solution["theta"](times),
                #parameters
            )
        except pybamm.SolverError:
            print('solver errored for params',parameters)
            return np.zeros_like(times)

    def n_parameters(self):
        return 6


if __name__ == '__main__':
    

    #FOR DIGIELCHM comparison
#=============================================================================
    files = [['C:/Users/natha/Desktop/Code/DigiElech/2023-06-06 Solution only/SC/SC_k0_1e-3_Ds_1e-5_Dp_1e-5.txt', 1e-3, 1e-5]]
    for i in files:
        print(i)
        
        volt = []
        curr = []
        row = []
        count = 0 
        f = open(i[0],'r')
        for row in f:
            count += 1
            row = row.split("\t")
            print(row)
            if row[0] == '':
                continue
            else:
                volt.append(float(row[0]))
                curr.append(float(row[1]))
#=============================================================================
            
            
        # pybamm.set_logging_level('INFO')
        model = SingleReactionSolution()
        
        React_Rate = 100
        Revers_Potent = 0
        Symm_fact = 0.5
        Uncomp_resist = 0
        Capacit = 0
        Volt_freq = 0
        x = np.array([React_Rate, Revers_Potent, Symm_fact, Uncomp_resist, Capacit, Volt_freq])
        
        Tmax = abs(0.5 + 0.5)/0.05 * 2
        Tdim = np.linspace(0, Tmax, 2**8)
        TnonDim = (96485.3328959 * 0.05 / (8.314459848*298)) * Tdim
        Tmax_nd = (96485.3328959 * 0.05 / (8.314459848*298)) * Tmax
        m = 2**8
        deltaT_nd = Tmax_nd / m

        t0 = time.perf_counter()
        current, voltage, theta = model.simulate(x, TnonDim)
        t1 = time.perf_counter()
        
        #faradaic current
        I_f = []
        I_f.append(0)
        
        for v in range(1, len(theta)):
            
            #generate faradaic current
            dOdt = (theta[v]-theta[v-1])/deltaT_nd
            I =  dOdt
            I_f.append(I)
                
        current = np.array(I_f)
        curr = np.array(curr)
                
        ##redimensionalizing here for now. Messy to do in main, move later
        I_d = current
        E_d = voltage * model._E_0
        #title = i[0]
        
        plt.cla()
        plt.plot(volt, curr/curr[0], color = 'g', label = 'Digielch')
        plt.plot(E_d, theta, color = 'b', label = 'Pybamm')
        #plt.title(title)
        plt.legend()
        plt.xlabel("Eapp [V]")
        plt.ylabel("current [A]")
        plt.show()
        #plt.savefig(i[0]+".png")

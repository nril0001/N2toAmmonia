## Basic Redox in Bulk Solution - diffusion controlled

import pybamm
import numpy as np

class CatalyticModel:
    def __init__(self,const_parameters,seioptions):

        #Const_parameters are the parameters to be passed into this model via main.py
        param = pybamm.ParameterValues(const_parameters)

        #Options relating to the models, like SEI
        self.options= self.calloptions()

        # Create dimensional fixed parameters
        # the "_d" indicates the dimensional form, while the lack of one/prescence of "_nd" means it's nondimensional
        #Adding separate D parameters for S and P
        CS_d = pybamm.Parameter("Far-field concentration of S(soln) [mol cm-3]")
        CP_d = pybamm.Parameter("Far-field concentration of P(soln) [mol cm-3]")
        DS_d = pybamm.Parameter("Diffusion Coefficient of S [cm2 s-1]")
        DP_d = pybamm.Parameter("Diffusion Coefficient of P [cm2 s-1]")
        F = pybamm.Parameter("Faraday Constant [C mol-1]")
        R = pybamm.Parameter("Gas constant [J K-1 mol-1]")
        a = pybamm.Parameter("Electrode Area [cm2]")
        T = pybamm.Parameter("Temperature [K]")
        
        E_start_d = pybamm.Parameter("Voltage start [V]")
        E_reverse_d = pybamm.Parameter("Voltage reverse [V]")
        v = pybamm.Parameter("Scan Rate [V s-1]")

        # Create dimensional input parameters
        E0_d = pybamm.InputParameter("Reversible Potential [V]")
        k0_d = pybamm.InputParameter("Redox Rate [s-1]")
        alpha = pybamm.InputParameter("Symmetry factor [non-dim]")
        Cdl_d = pybamm.InputParameter("Capacitance [F]")
        Ru_d = pybamm.InputParameter("Uncompensated Resistance [Ohm]")

        # Create scaling factors   for non-dimensionalisation
        Dmax = DS_d + DP_d
        E_0 = (R * T) / F #units are V
        T_0 = E_0 / v #units are seconds
        X_0 = pybamm.sqrt((F*v)/(R*T*Dmax)) #units are cm
        #I_0 = (DS_d * F * a * CS_d)/ X_0  #units are A
        I_0 = F*a*CS_d*pybamm.sqrt(Dmax)*pybamm.sqrt((F*v)/(R*T))
        K_0 = pybamm.sqrt((R*T*Dmax)/(F*v))/DS_d
        Cdl_0 = (a * E_0)/(I_0 * T_0) # V/A s
        Ru_0 = I_0 / E_0 #units are Amps/V
        
        # Non-dimensionalise parameters
        E0 = E0_d / E_0 #no units
        k0 = k0_d * K_0 #no units
        Cdl = Cdl_d * Cdl_0 #no units
        Ru = Ru_d * Ru_0 #no units
        
        #Diffusion coefficients
        d_S = DS_d/Dmax
        d_P = DP_d/Dmax
        
        #Concentrations
        Ctot = CS_d + CP_d
        cs_nd = (CS_d/Ctot)
        cp_nd = (CP_d/Ctot)

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
        Edc_backwards = pybamm.t - 2 * t_reverse
        Eapp = E_start + \
            (pybamm.t <= t_reverse) * Edc_forward + \
            (pybamm.t > t_reverse) * Edc_backwards                   

        # create PyBaMM model object
        model = pybamm.BaseBatteryModel(options=seioptions)

        # Create state variables for model
        #sc is surface concentration
        c_s = pybamm.Variable("S(soln) [non-dim]", domain="solution")
        c_p = pybamm.Variable("P(soln) [non-dim]", domain="solution")
        i = pybamm.Variable("Current [non-dim]")

        # Effective potential
        Eeff = Eapp - i * Ru #no units
        

        # defining boundary values for S and P
        c_at_electrode_s = pybamm.BoundaryValue(c_s, "left")
        c_at_electrode_p = pybamm.BoundaryValue(c_p, "left")

        # Faradaic current (Butler Volmer)
        #i_f = pybamm.BoundaryGradient(c_s, "left")
        butler_volmer = k0 * ((c_at_electrode_p) * pybamm.exp((1 - alpha) * (Eeff - E0))
                            - ((c_at_electrode_s) * pybamm.exp(-alpha * (Eeff - E0))))  

        #"left" indicates environment directly on electrode surface; x = 0
        #"right" indicates environment between diffusion layer and bulk solution; x = xmax 

        # PDEs - left hand side is assumed to be time derivative of the PDE
        #multiplying by the other coefficient...? idk why but it works
        model.rhs = {
            c_s: d_P * pybamm.div(pybamm.grad(c_s)),
            c_p: d_S * pybamm.div(pybamm.grad(c_p)),
            i: butler_volmer - i,
        }
        
        # Setting boundary and initial conditions
        model.boundary_conditions = {
            c_s: {
                "right": (cs_nd, "Dirichlet"),
                "left": (-butler_volmer, "Neumann"),                    
            },

            c_p: {
                "right": (cp_nd, "Dirichlet"),
                "left": (butler_volmer, "Neumann"),                 
            } 
        }

        model.initial_conditions = {
            c_s: cs_nd,
            c_p: cp_nd,
            i: Cdl * (1.0),
        }

        # set spatial variables and domain geometry
        x = pybamm.SpatialVariable('x', domain="solution")
        x_max = 6*pybamm.sqrt(d_S*Tmax_nd)
        model.geometry = pybamm.Geometry({
            "solution": {
                    x: {
                        "min": pybamm.Scalar(0),
                        "max": x_max
                    }
            }
        })

        # Controls how many space steps are taken (higher number, higher accuracy)
        model.var_pts = {
            x: 100
        }

        # Using Finite Volume discretisation on an expanding 1D grid
        model.submesh_types = {
            "solution": pybamm.MeshGenerator(
                pybamm.Exponential1DSubMesh, {'side': 'left'}
            )
        }
        model.spatial_methods = {
            "solution": pybamm.FiniteVolume()
        }

        # model variables
        model.variables = {
            "Current [non-dim]": i,
            "Applied Voltage [non-dim]": Eapp,
            "S(soln) at electrode [non-dim]": c_at_electrode_s,
            "P(soln) at electrode [non-dim]": c_at_electrode_p,
        }
        
        # Set model parameters
        param.process_model(model)
        geometry = model.geometry
        param.process_geometry(geometry)
        

        # Create mesh and discretise model
        mesh = pybamm.Mesh(
            geometry, model.submesh_types, model.var_pts)
        disc = pybamm.Discretisation(mesh, model.spatial_methods)
        disc.process_model(model)

        # Create solver
        model.convert_to_format = 'python'
        solver = pybamm.ScipySolver(method='Radau', rtol=1e-6, atol=1e-6)
        # model.convert_to_format = 'casadi'
        # solver = pybamm.CasadiSolver(mode='safe', rtol=1e-9, atol=1e-9, root_method='casadi')
        
    
        # Store discretised model and solver
        self._model = model
        self._param = param
        self._solver = solver

        #Store processed symbols/dimensionalization factors
        self._E_0 = param.process_symbol(E_0).evaluate()
        self._I_0 = param.process_symbol(I_0).evaluate()
        self._T_0 = param.process_symbol(T_0).evaluate()
        self._CS_d = param.process_symbol(CS_d).evaluate()
        self._CP_d = param.process_symbol(CP_d).evaluate()
        self._deltaT_nd = param.process_symbol(deltaT_nd).evaluate()

        # store time scale related things
        self._Tmax_nd = param.process_symbol(Tmax_nd).evaluate()
        self._m = m
        
        print("Catalytic Model 04 initialized successfully.")

    def simulate(self, parameters):
        #####DEBUGGING#####
        #pybamm.set_logging_level("DEBUG")

        #7 May 23: method to pull times from init
        times_nd = np.linspace(0, self._Tmax_nd, int(self._m))
        print(f"Number of timesteps: " + str(self._m))
        try:
            solution = self._solver.solve(self._model, times_nd, inputs=parameters)
            c_S = solution["S(soln) at electrode [non-dim]"](times_nd)
            c_P = solution["P(soln) at electrode [non-dim]"](times_nd)
            E = solution["Applied Voltage [non-dim]"](times_nd)
            
            #faradaic current
            I_f = []
            I_f.append(0)
            
            #capacitative current
            I_cdl = []
            I_cdl.append(0)
            
            for v in range(1, self._m):
                
                #generate faradaic current
                dOdt = (c_S[v]-c_S[v-1])/self._deltaT_nd
                I =  dOdt
                I_f.append(I)
                
            current = []
            for v in range(0, len(I_f)):
                if v < len(I_f)/2:
                    current.append(I_f[v])
                else:
                    current.append(I_f[v])
                    
            current = np.array(current)
            
        except pybamm.SolverError as e:
            print(e)
            solution = np.zeros_like(times_nd)
        return (
            current, E,
            c_S, c_P,
            times_nd,
        )

#TODO: Make a redimensionalise function
    # nd_sol are nondimensional solutions
    def redimension(self, nd_sol):
        #need to get factors in here to multiply
        #I_nd, E_nd, O_nd, S_nd, P_nd = cmodel.simulate(input_parameters,times)
        #
        return

    #FIXME: make a call function for the options here
    def calloptions(seioptions):
       return
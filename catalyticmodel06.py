## Basic Redox in Bulk Solution - diffusion controlled

import pybamm
import numpy as np

class CatalyticModel:
    def __init__(self,const_parameters,seioptions, atoler, rtoler, t_steps, x_steps):

        #Const_parameters are the parameters to be passed into this model via main.py
        self.const_parameters = const_parameters
        self.seioptions = seioptions
        self.atoler = atoler
        self.rtoler = rtoler
        self.t_steps = t_steps
        self.x_steps = x_steps

        #Options relating to the models, like SEI
        self.options= self.calloptions()

        # Create dimensional fixed parameters
        # the "_d" indicates the dimensional form, while the lack of one/prescence of "_nd" means it's nondimensional
        #Adding separate D parameters for S and P
        self.CS_d = pybamm.Parameter("Far-field concentration of S(soln) [mol cm-3]")
        self.CP_d = pybamm.Parameter("Far-field concentration of P(soln) [mol cm-3]")
        self.DS_d = pybamm.Parameter("Diffusion Coefficient of S [cm2 s-1]")
        self.DP_d = pybamm.Parameter("Diffusion Coefficient of P [cm2 s-1]")
        self.Delta = pybamm.Parameter("Diffusion Layer Thickness [cm]")
        self.F = pybamm.Parameter("Faraday Constant [C mol-1]")
        self.R = pybamm.Parameter("Gas constant [J K-1 mol-1]")
        self.a = pybamm.Parameter("Electrode Area [cm2]")
        self.r = pybamm.Parameter("Electrode Radius [cm]")
        self.T = pybamm.Parameter("Temperature [K]")
        
        self.E_start_d = pybamm.Parameter("Voltage start [V]")
        self.E_reverse_d = pybamm.Parameter("Voltage reverse [V]")
        self.v = pybamm.Parameter("Scan Rate [V s-1]")

        # Create dimensional input parameters
        self.E0_d = pybamm.InputParameter("Reversible Potential [V]")
        self.k0_d = pybamm.InputParameter("Redox Rate [cm s-1]")
        self.kads_d = pybamm.InputParameter("Adsorption Rate [mol-1 cm3 s-1]")
        self.kdes_d = pybamm.InputParameter("Desorption Rate [mol-1 cm3 s-1]")
        self.k1_d = pybamm.InputParameter("Catalytic Rate [s-1]")
        self.alpha = pybamm.InputParameter("Symmetry factor [non-dim]")
        self.Cdl_d = pybamm.InputParameter("Capacitance [F]")
        self.Ru_d = pybamm.InputParameter("Uncompensated Resistance [Ohm]")

        # Create scaling factors   for non-dimensionalisation
        self.T_0 = self.DS_d/self.r**2 #units are s-1
        self.X_0 = 1/self.r #units are cm-1
        self.K_0 = self.r/self.DS_d #units are s/cm, heterogenous
        self.KADS_0 = 1
        self.KDES = 1
        self.K_1 = self.r**2/self.DS_d #units are s-1, homogenous
        self.C_0 = 1/self.CS_d # concentration
        self.D_0 = 1/self.DS_d # diffusion
        self.E_0 = self.F / (self.R * self.T)  #units are V-1
        self.V_0 = (self.r**2/self.DS_d)*(self.F / (self.R * self.T)) #scan rate
        self.I_0 = 1/(np.pi*self.r*self.F*self.DS_d*self.CS_d)# current A-1
        
        self.Cdl_0 = (self.I_0 * self.T_0)/(self.E_0) # V/A s
        self.Ru_0 = self.I_0/self.E_0 #units are Amps/V
        
        # Non-dimensionalise parameters
        self.E0 = self.E0_d * self.E_0 #no units
        self.k0 = self.k0_d * self.K_0 #no units
        self.k1 = self.k1_d * self.K_1 # no units
        self.kads = self.kads_d * self.KADS_0 # no units
        self.kdes = self.kdes_d * self.KDES_0 # no units
        self.Cdl = self.Cdl_d * self.Cdl_0 #no units
        self.Ru = self.Ru_d * self.Ru_0 #no units
        
        #Diffusion coefficients
        self.d_S = self.DS_d * self.D_0
        self.d_P = self.DP_d * self.D_0
        self.d_max = pybamm.maximum(self.d_S, self.d_P)
        
        #Concentrations
        self.cs_nd = self.CS_d * self.C_0
        self.cp_nd = self.CP_d * self.C_0

        self.E_start = self.E_start_d * self.E_0
        self.E_reverse = self.E_reverse_d * self.E_0
        
        #creating time scale and non-dimensionalizing
        self.Tmax_nd = (abs(self.E_start_d - self.E_reverse_d)/self.v)*self.T_0
        
        #length of time step, nondimensional
        self.deltaT_nd = self.Tmax_nd / self.t_steps
        
        #max length of diffusion layer
        self.x_max = 6 * pybamm.sqrt(self.d_max * self.Tmax_nd) 
        # self.x_max = self.Delta
        
        #scan rate
        self.V = self.V_0 * self.v
        
    def model(self, sweep="forward", Cs =None, Cp =None, Ee=None, Ea=None, Cs_a=None, Cp_a=None):
        
        param = pybamm.ParameterValues(self.const_parameters)
        
        # create PyBaMM model object
        model = pybamm.BaseBatteryModel(options=self.seioptions)

        # Create state variables for model
        #sc is surface concentration
        c_s = pybamm.Variable("S(soln) [non-dim]", domain="solution")
        c_p = pybamm.Variable("P(soln) [non-dim]", domain="solution")
        c_s_a = pybamm.Variable("S(ads) [non-dim]")
        c_p_a = pybamm.Variable("P(ads) [non-dim]")
        Eeff = pybamm.Variable("Effective Voltage [non-dim]")
        
        if sweep == "forward":
            Eapp = self.E_start - self.V * pybamm.t
            Cs = self.cs_nd
            Cp = self.cp_nd
            Cs_a = 0
            Cp_a = 0
            Ee = self.E_start
            Ea = self.E_start
        elif sweep == "backward":
            Eapp = self.E_reverse + self.V * pybamm.t
            Cs = pybamm.Array(Cs, domain="solution")
            Cp = pybamm.Array(Cp, domain="solution")
            Cs_a = Cs_a
            Cp_a = Cp_a
            Ee = Ee 
            Ea = Ea

        # defining boundary values for S and P
        c_at_electrode_s = pybamm.BoundaryValue(c_s, "left")
        c_at_electrode_p = pybamm.BoundaryValue(c_p, "left")

        # Faradaic current (Butler Volmer)
        bv1 = self.k0 * ((c_at_electrode_p) * pybamm.exp((1 - self.alpha) * (Eeff - self.E0))
                            - ((c_at_electrode_s) * pybamm.exp(-self.alpha * (Eeff - self.E0))))
        
        bv2 = self.k0 * ((c_at_electrode_p) * pybamm.exp((1 - self.alpha) * (Eeff - self.E0))
                            - ((c_at_electrode_s) * pybamm.exp(-self.alpha * (Eeff - self.E0))))
        
        dOdt = bv1 + bv2
        
        i_f = dOdt
        
        i_cap = self.Cdl
        
        # csa = (-self.kasd*c_s_a) 
        # cpa =
        
        
        if sweep == "forward":
            i = i_f - i_cap
        elif sweep == "backward":
            i = i_f + i_cap 

        # PDEs - left hand side is assumed to be time derivative of the PDE
        #dividing by their own coefficients
        model.rhs = {
            c_s_a: 1,
            c_p_a: 1,
            c_s: (pybamm.div(pybamm.grad(c_s)) * self.d_S),
            c_p: (pybamm.div(pybamm.grad(c_p)) * self.d_P),
        }
        
        # algebraic equations (none)
        model.algebraic = {
            Eeff: Eapp - self.Ru*i - Eeff,
        }
        
        # Setting boundary and initial conditions
        #"left" indicates environment directly on electrode surface; x = 0
        #"right" indicates environment between diffusion layer and bulk solution; x = xmax
        model.boundary_conditions = {
            c_s: {
                "right": (self.cs_nd, "Dirichlet"),
                "left": ((-bv1/self.d_S), "Neumann"),                    
            },

            c_p: {
                "right": (self.cp_nd, "Dirichlet"),
                "left": ((bv1/self.d_P), "Neumann"),                 
            } 
        }
        
        model.initial_conditions = {
            c_s_a: Cs_a,
            c_p_a: Cp_a,
            c_s: Cs,
            c_p: Cp,
            Eeff: Ee,
        }

        # set spatial variables and domain geometry
        x = pybamm.SpatialVariable('x', domain="solution")
        model.geometry = pybamm.Geometry({
            "solution": {
                    x: {
                        "min": pybamm.Scalar(0),
                        "max": self.x_max
                    }
            }
        })

        # Controls how many space steps are taken (higher number, higher accuracy)
        model.var_pts = {
            x: self.x_steps
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
            "Effective Voltage [non-dim]": Eeff,
            "S(soln) at electrode [non-dim]": c_at_electrode_s,
            "P(soln) at electrode [non-dim]": c_at_electrode_p,
            "S(soln) [non-dim]": c_s,
            "P(soln) [non-dim]": c_p,
            
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
        solver = pybamm.ScikitsDaeSolver(atol=self.atoler, rtol=self.rtoler)
    
        # Store discretised model and solver
        self._model = model
        self._param = param
        self._solver = solver

        #Store processed symbols/dimensionalization factors
        self._E_0 = param.process_symbol(self.E_0).evaluate()
        self._I_0 = param.process_symbol(self.I_0).evaluate()
        self._T_0 = param.process_symbol(self.T_0).evaluate()
        self._CS_d = param.process_symbol(self.CS_d).evaluate()
        self._CP_d = param.process_symbol(self.CP_d).evaluate()
        self._deltaT_nd = param.process_symbol(self.deltaT_nd).evaluate()

        # store time scale related things
        self._Tmax_nd = param.process_symbol(self.Tmax_nd).evaluate()
        self._m = self.t_steps
        self._x = self.x_steps
        
        print("Catalytic Model 05 initialized successfully.")

    def simulate(self, parameters):
        #####DEBUGGING#####
        # pybamm.set_logging_level("DEBUG")

        #7 May 23: method to pull times from init
        self.model("forward")
        times_nd = np.linspace(0, self._Tmax_nd, int(self._m)//2)
        times = np.linspace(0, self._Tmax_nd*2, int(self._m))
        print("Number of timesteps: " + str(self._m))
        print("Number of spacesteps: " + str(self._x))
        try:
            solution = self._solver.solve(self._model, times_nd, inputs=parameters)
            c_S = solution["S(soln) [non-dim]"](times_nd)
            c_P = solution["P(soln) [non-dim]"](times_nd)
            cS = solution["S(soln) at electrode [non-dim]"](times_nd)
            cP = solution["P(soln) at electrode [non-dim]"](times_nd)
            E = solution["Applied Voltage [non-dim]"](times_nd)
            Ee = solution["Effective Voltage [non-dim]"](times_nd)
            current = solution["Current [non-dim]"](times_nd)
            time = solution["time [non-dim]"](times_nd)
            
     
            self.model("backward", c_S[1:-1, -1], c_P[1:-1, -1], Ee[-1], E[-1])
            solution = self._solver.solve(self._model, times_nd, inputs=parameters)
            c_S = np.concatenate((c_S, solution["S(soln) [non-dim]"](times_nd)))
            c_P = np.concatenate((c_P, solution["P(soln) [non-dim]"](times_nd)))
            cS = np.concatenate((cS, solution["S(soln) at electrode [non-dim]"](times_nd)))
            cP = np.concatenate((cP, solution["P(soln) at electrode [non-dim]"](times_nd)))
            E = np.concatenate((E, solution["Applied Voltage [non-dim]"](times_nd)))
            Ee = np.concatenate((Ee, solution["Effective Voltage [non-dim]"](times_nd)))
            current = np.concatenate((current, solution["Current [non-dim]"](times_nd)))
            
        except pybamm.SolverError as e:
            print(e)
            solution = np.zeros_like(times_nd)
            
        return (
            current, E,
            cS, cP,
            times
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
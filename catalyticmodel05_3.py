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
        self.F = pybamm.Parameter("Faraday Constant [C mol-1]")
        self.R = pybamm.Parameter("Gas constant [J K-1 mol-1]")
        self.a = pybamm.Parameter("Electrode Area [cm2]")
        self.T = pybamm.Parameter("Temperature [K]")
        
        self.E_start_d = pybamm.Parameter("Voltage start [V]")
        self.E_reverse_d = pybamm.Parameter("Voltage reverse [V]")
        self.v = pybamm.Parameter("Scan Rate [V s-1]")

        # Create dimensional input parameters
        self.E0_d = pybamm.InputParameter("Reversible Potential [V]")
        self.k0_d = pybamm.InputParameter("Redox Rate [s-1]")
        self.alpha = pybamm.InputParameter("Symmetry factor [non-dim]")
        self.Cdl_d = pybamm.InputParameter("Capacitance [F]")
        self.Ru_d = pybamm.InputParameter("Uncompensated Resistance [Ohm]")

        # Create scaling factors   for non-dimensionalisation
        self.E_0 = (self.R * self.T) / self.F #units are V
        self.T_0 = self.E_0 / self.v #units are seconds
        self.X_0 = pybamm.sqrt((self.F*self.v)/(self.R*self.T*self.DS_d)) #units are cm
        #I_0 = (DS_d * F * a * CS_d)/ X_0  #units are A
        self.I_0 = self.F*self.a*self.CS_d*pybamm.sqrt(self.DS_d)*pybamm.sqrt((self.F*self.v)/(self.R*self.T))
        self.K_0 = (pybamm.sqrt(self.R * self.T * self.DS_d/self.F * self.v))/self.DS_d #units are s
        #K_0 = a / Dmax
        self.Cdl_0 = (self.a * self.E_0)/(self.I_0 * self.T_0) # V/A s
        self.Ru_0 = self.I_0 / self.E_0 #units are Amps/V
        
        # Non-dimensionalise parameters
        self.E0 = self.E0_d / self.E_0 #no units
        self.k0 = self.k0_d * self.K_0 #no units
        self.Cdl = self.Cdl_d * self.Cdl_0 #no units
        self.Ru = self.Ru_d * self.Ru_0 #no units
        
        #Diffusion coefficients
        self.d_S = self.DS_d/self.DS_d
        self.d_P = self.DP_d/self.DS_d
        self.d_max = pybamm.maximum(self.d_S, self.d_P)
        
        #Concentrations
        self.Ctot = self.CS_d + self.CP_d
        self.cs_nd = (self.CS_d/self.Ctot)
        self.cp_nd = (self.CP_d/self.Ctot)

        self.E_start = self.E_start_d / self.E_0
        self.E_reverse = self.E_reverse_d / self.E_0
        
        #creating time scale and non-dimensionalizing
        self.Tmax_nd = ((abs(self.E_start_d - self.E_reverse_d) / self.v * 2)/ self.T_0)
        
        #length of time step, nondimensional
        self.deltaT_nd = self.Tmax_nd / self.t_steps
        
        #max length of diffusion layer
        self.x_max = 6 * pybamm.sqrt(self.d_max * self.Tmax_nd)
        
    def model(self, sweep="forward", Cs =None, Cp =None, Ee=None, Ea=None):
        
        param = pybamm.ParameterValues(self.const_parameters)
        
        if sweep == "forward":
            Edc = -1
            i_cap = -self.Cdl * self.deltaT_nd
            Cs = self.cs_nd
            Cp = self.cp_nd
            Ee = self.E_start
            Ea = self.E_start
        elif sweep == "backward":
            Edc = 1
            i_cap = self.Cdl * self.deltaT_nd
            Cs = pybamm.Array(Cs, domain="solution")
            Cp = pybamm.Array(Cp, domain="solution")
            Ee = Ee 
            Ea = Ea
        
        # create PyBaMM model object
        model = pybamm.BaseBatteryModel(options=self.seioptions)

        # Create state variables for model
        #sc is surface concentration
        c_s = pybamm.Variable("S(soln) [non-dim]", domain="solution")
        c_p = pybamm.Variable("P(soln) [non-dim]", domain="solution")
        Eeff = pybamm.Variable("Effective Voltage [non-dim]")
        Eapp = pybamm.Variable("Applied Voltage [non-dim]")

        # defining boundary values for S and P
        c_at_electrode_s = pybamm.BoundaryValue(c_s, "left")
        c_at_electrode_p = pybamm.BoundaryValue(c_p, "left")

        # Faradaic current (Butler Volmer)
        butler_volmer = self.k0 * ((c_at_electrode_p) * pybamm.exp((1 - self.alpha) * (Eeff - self.E0))
                            - ((c_at_electrode_s) * pybamm.exp(-self.alpha * (Eeff - self.E0))))  
        
        dOdt = butler_volmer
        
        i_f = dOdt
        
        i = i_f + i_cap

        #"left" indicates environment directly on electrode surface; x = 0
        #"right" indicates environment between diffusion layer and bulk solution; x = xmax 

        # PDEs - left hand side is assumed to be time derivative of the PDE
        #dividing by their own coefficients
        model.rhs = {
            Eapp: Edc,
            c_s: pybamm.div(pybamm.grad(c_s)) * self.d_S,
            c_p: pybamm.div(pybamm.grad(c_p)) * self.d_P,
        }
        
        # algebraic equations (none)
        model.algebraic = {
            Eeff: Eapp - self.Ru*i - Eeff,
        }
        
        # Setting boundary and initial conditions
        model.boundary_conditions = {
            c_s: {
                "right": (self.cs_nd, "Dirichlet"),
                "left": ((-butler_volmer/self.d_S), "Neumann"),                    
            },

            c_p: {
                "right": (self.cp_nd, "Dirichlet"),
                "left": ((butler_volmer/self.d_P), "Neumann"),                 
            } 
        }
        
        model.initial_conditions = {
            c_s: Cs,
            c_p: Cp,
            Eeff: Ee,
            Eapp: Ea,
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
            "P(soln) [non-dim]": c_p
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
        # model.convert_to_format = 'python'
        solver = pybamm.ScikitsDaeSolver(atol=self.atoler, rtol=self.rtoler)
        # model.convert_to_format = 'casadi'
        # solver = pybamm.IDAKLUSolver(atol=atoler, rtol=rtoler)
        # solver = pybamm.CasadiSolver(mode='fast', rtol=rtoler, atol=atoler, root_method='casadi')
    
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
        times_nd = np.linspace(0, self._Tmax_nd/2, int(self._m)//2)
        print("Number of timesteps: " + str(self._m))
        print("Number of spacesteps: " + str(self._x))
        try:
            solution = self._solver.solve(self._model, times_nd, inputs=parameters)
            c_S_f = solution["S(soln) [non-dim]"](times_nd)
            c_P_f = solution["P(soln) [non-dim]"](times_nd)
            E_f = solution["Applied Voltage [non-dim]"](times_nd)
            Ee_f = solution["Effective Voltage [non-dim]"](times_nd)
            current_f = solution["Current [non-dim]"](times_nd)
            
            #removes the boundary conditions from the final output array (theoretically)
            #don't know if the boundary conditions will be overwritten if an array of len 102 is used
            cS = c_S_f[1:-1, -1]
            cP = c_P_f[1:-1, -1]
            yes =c_S_f[:, -1]
            no = c_P_f[:, -1]
            maybe = Ee_f[-1]
            definitely = E_f[-1]
     
            self.model("backward", c_S_f[1:-1, -1], c_P_f[1:-1, -1], Ee_f[-1], E_f[-1])
            solution = self._solver.solve(self._model, times_nd, inputs=parameters)
            c_S_b = solution["S(soln) [non-dim]"](times_nd)
            c_P_b = solution["P(soln) [non-dim]"](times_nd)
            E_b = solution["Applied Voltage [non-dim]"](times_nd)
            Ee_b = solution["Effective Voltage [non-dim]"](times_nd)
            current_b = solution["Current [non-dim]"](times_nd)
            
            # print(current_f[-1], current_b[1])
            
            c_S = np.concatenate((c_S_f, c_S_b))
            c_P = np.concatenate((c_P_f, c_P_b))
            E = np.concatenate((E_f, E_b))
            Ee = np.concatenate((Ee_f, Ee_b))
            current = np.concatenate((current_f, current_b))
            
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
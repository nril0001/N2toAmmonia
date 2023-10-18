## Adsorption concerted electron transfer 
## Single layer lithium deposition and Li3N reaction only
## Li+ + S* + e- <-> Li*
## 3Li* + N <-> Li3N* (simplified to Li* + X <-> LiX*)

import pybamm
import numpy as np

class CatalyticModel:
    def __init__(self,const_parameters,seioptions, atoler, rtoler, t_steps, x_steps, solver):

        #Const_parameters are the parameters to be passed into this model via main.py
        self.const_parameters = const_parameters
        self.seioptions = seioptions
        self.atoler = atoler
        self.rtoler = rtoler
        self.t_steps = t_steps
        self.x_steps = x_steps
        self._solver_ = solver

        #Options relating to the models, like SEI
        self.options= self.calloptions()

        # Create dimensional fixed parameters
        # the "_d" indicates the dimensional form, while the lack of one/prescence of "_nd" means it's nondimensional
        #Adding separate D parameters for S and P
        self.CS_d = pybamm.Parameter("Far-field concentration of S(soln) [mol cm-3]")
        self.CY_d = pybamm.Parameter("Far-field concentration of Y(soln) [mol cm-3]")
        self.C0 = pybamm.Parameter("Standard Unity Concentration [mol cm-3]")
        
        self.SCP_d = pybamm.Parameter("Surface coverage of P [mol cm-2]")
        self.SCX_d = pybamm.Parameter("Surface coverage of X [mol cm-2]")
        self.SCZ_d = pybamm.Parameter("Surface coverage of Z [mol cm-2]")
        
        self.DS_d = pybamm.Parameter("Diffusion Coefficient of S [cm2 s-1]")
        self.DY_d = pybamm.Parameter("Diffusion Coefficient of Y [cm2 s-1]")
        
        self.Delta = pybamm.Parameter("Diffusion Layer Thickness [cm]")
        self.Gamma = pybamm.Parameter("Electrode Coverage [mol cm-2]")
        self.a = pybamm.Parameter("Electrode Area [cm2]")
        self.r = pybamm.Parameter("Electrode Radius [cm]")
        
        self.F = pybamm.Parameter("Faraday Constant [C mol-1]")
        self.R = pybamm.Parameter("Gas constant [J K-1 mol-1]")
        self.T = pybamm.Parameter("Temperature [K]")
        
        self.Cdl_d = pybamm.Parameter("Capacitance [F]")
        self.Ru_d = pybamm.Parameter("Uncompensated Resistance [Ohm]")
        
        self.E_start_d = pybamm.Parameter("Voltage start [V]")
        self.E_reverse_d = pybamm.Parameter("Voltage reverse [V]")
        self.v = pybamm.Parameter("Scan Rate [V s-1]")

        # Create dimensional input parameters"
        self.E01_d = pybamm.InputParameter("Reversible Potential 1 [V]")
        self.k0_d = pybamm.InputParameter("Electrosorption Rate [mol-1 cm3 s-1]")
        self.k1_d = pybamm.InputParameter("Catalytic Rate For [mol-1 cm3 s-1]")
        self.k1b_d = pybamm.InputParameter("Catalytic Rate Back [s-1]")
        self.alpha = pybamm.InputParameter("Symmetry factor [non-dim]")
        self.G = pybamm.InputParameter("G")
        self.G_ = pybamm.InputParameter("G'")

        # Create scaling factors for non-dimensionalisation
        self.T_0 = self.DS_d/self.r**2 # time, units in s-1
        self.X_0 = 1/self.r # distance, units in cm-1
        self.C_0 = 1/self.CS_d # concentration, units in mol cm-3 
        self.G_0 = 1/self.Gamma # surface concentration units in mol-1 cm2
        self.K_0ads = (self.r**2)/self.DS_d # redox rate constant (ads), units in s-1
        self.K_0sol = (self.r)/self.DS_d # redox rate constant (sol), units in mol cm s
        self.K0ads_0 = (self.Gamma*self.r)/self.DS_d # adsorption rate constant, units in mol-3 cm s
        self.Kdes_0 = (self.Gamma*self.r)/(self.DS_d*self.CS_d) # adsorption rate constant, units in s
        self.D_0 = 1/self.DS_d # diffusion, units in s cm-2
        self.V_0 = (self.r**2/self.DS_d)*(self.F / (self.R * self.T)) # scan rate, units in s V-1
        self.E_0 = self.F / (self.R * self.T) # potential, units in V-1
        self.B_0 = (self.r*self.CS_d)/self.Gamma #Saturation coefficient
        self.I_0 = 1/(np.pi*self.r*self.F*self.DS_d*self.CS_d) # current, units in A-1
        
        self.Cdl_0 = (self.E_0)/(self.I_0 * self.T_0) # capacitance, units in V s A-1
        self.Ru_0 = self.E_0/self.I_0 # resistance, units in A V-1
        
        # Non-dimensionalise parameters
        #Rate constants
        self.k0 = self.k0_d * self.K0ads_0 #no units
        self.k1 = self.k1_d * self.K0ads_0
        self.k1b = self.k1b_d * self.Kdes_0
        
        #Diffusion coefficients
        self.d_S = self.DS_d * self.D_0 #no units
        self.d_Y = self.DY_d * self.D_0 #no units
        self.d_max = pybamm.maximum(self.d_S, self.d_Y) #no units
        
        #Concentrations
        self.c0 = self.C0 * self.C_0 #no units
        self.cs_nd = self.CS_d * self.C_0 #no units
        self.cy_nd = self.CY_d * self.C_0 #no units
        
        #Surface concentrations
        self.G_max = self.Gamma*self.G_0 #no units
        self.scp_nd = self.SCP_d * self.G_0 # no units
        self.scx_nd = self.SCX_d * self.G_0 # no units
        self.scz_nd = self.SCZ_d * self.G_0 # no units

        #Potential
        self.E01 = self.E01_d * self.E_0 #no units
        self.E_start = self.E_start_d * self.E_0 #no units
        self.E_reverse = self.E_reverse_d * self.E_0 #no units
        
        #Scan Rate
        self.V = self.V_0 * self.v #no units
        
        #Time
        self.Tmax_nd = (abs(self.E_start_d - self.E_reverse_d)/self.v)*self.T_0 #no units
        self.deltaT_nd = self.Tmax_nd / self.t_steps #no units
        
        #Distance
        # self.x_max = 6 * pybamm.sqrt(self.d_max * self.Tmax_nd) #no units
        self.x_max = self.Delta * self.X_0
        
        #Resistance and Capacitance
        self.Cdl = self.Cdl_d * self.Cdl_0 #no units
        self.Ru = self.Ru_d * self.Ru_0 #no units
    
    #define concentration of reactant electrons
    def BV_red(self, E, E0):
        return pybamm.exp((-self.alpha) * (E - E0))
    
    #define concentration of product electrons
    def BV_ox(self, E, E0):
        return pybamm.exp((1 - self.alpha) * (E - E0))
     
    def model(self, sweep="forward", CS=None, CY=None, SCP=None, SCX=None, SCZ=None,  Ee=None):
        
        # create PyBaMM model object
        param = pybamm.ParameterValues(self.const_parameters)
        model = pybamm.BaseBatteryModel(options=self.seioptions)

        # Create state variables for model
        #sc is surface concentration
        c_s = pybamm.Variable("S(soln) [non-dim]", domain="solution")
        c_y = pybamm.Variable("Y(soln) [non-dim]", domain="solution")
        sc_p = pybamm.Variable("P(ads) [non-dim]")
        sc_x = pybamm.Variable("X(ads) [non-dim]")
        sc_z = pybamm.Variable("Z(ads) [non-dim]")
        Eeff = pybamm.Variable("Effective Voltage [non-dim]")
        
        #setting initial conditions
        if sweep == "forward":
            Eapp = self.E_start - self.V * pybamm.t
            Ee = self.E_start
            Cs = self.cs_nd
            Cy = self.cy_nd
            SCp = self.scp_nd
            SCx = self.G_max
            SCz = self.scz_nd
        elif sweep == "backward":
            Eapp = self.E_reverse + self.V * pybamm.t
            Ee = Ee
            Cs = pybamm.Array(CS, domain="solution")
            Cy = pybamm.Array(CY, domain="solution")
            SCp = SCP
            SCx = SCX
            SCz = SCZ

        # defining boundary values for S
        c_at_electrode_s = pybamm.BoundaryValue(c_s, "left")
        c_at_electrode_y = pybamm.BoundaryValue(c_y, "left")
        dOdx = pybamm.BoundaryGradient(c_s, "left")

        BV_red1 = self.BV_red(Eeff, self.E01)
        BV_ox1 = self.BV_ox(Eeff, self.E01)

        # Faradaic current (Butler Volmer)
        BV1 = self.k0 * ((c_at_electrode_s *  BV_red1 *pybamm.exp((-self.G_)*(sc_p))) 
                        - ((sc_p) * self.c0 * BV_ox1*pybamm.exp((self.G-self.G_)*(sc_p)))) 
        
        #Catalytic Reaction - 6Li + N2 
        cat_for = (self.k1 * c_at_electrode_y * sc_p) 
        cat_back = (self.k1b * sc_z)
        
        #time derivatives
        dSdt = (pybamm.div(pybamm.grad(c_s)) * self.d_S) #Lithium
        dYdt = (pybamm.div(pybamm.grad(c_y)) * self.d_Y) #Nitrogen
        
        dsPdt = (BV1 - cat_for + cat_back)*self.B_0 #Adsorbed Lithium
        dsZdt = (cat_for - cat_back)*self.B_0 # Lithium Nitride
        dsXdt = 1 # Active Sites
        
        #space derivatives
        dSdx = (BV1)/self.d_S
        dYdx = (cat_for - cat_back)/self.d_Y
        
        #current
        i_f = -dOdx
        i_cap = self.Cdl * self.V
        if sweep == "forward":
            i = i_f - i_cap
        elif sweep == "backward":
            i = i_f + i_cap

        #"left" indicates environment directly on electrode surface; x = 0
        #"right" indicates environment between diffusion layer and bulk solution; x = xmax 

        # PDEs - left hand side is assumed to be time derivative of the PDE
        #dividing by their own coefficients
        model.rhs = {
            c_s: dSdt,
            c_y: dYdt,
            sc_p: dsPdt,
            sc_x: dsXdt,
            sc_z: dsZdt,
        }
        
        # algebraic equations (none)
        model.algebraic = {
            Eeff: Eapp - self.Ru*i - Eeff,
        }
        
        # Setting boundary and initial conditions
        model.boundary_conditions = {
            c_s: {"right": (self.cs_nd, "Dirichlet"),
                  "left": ((dSdx), "Neumann"),},
            c_y: {"right": (self.cy_nd, "Dirichlet"),
                  "left": ((dYdx), "Neumann"),},
        }
        
        model.initial_conditions = {
            Eeff: Ee,
            c_s: Cs,
            c_y: Cy,
            sc_p: SCp,
            sc_x: SCx,
            sc_z: SCz,
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
            "S(soln) [non-dim]": c_s,
            "Y(soln) [non-dim]": c_y,
            "P(ads) [non-dim]": sc_p,
            "X(ads) [non-dim]": sc_x,
            "Z(ads) [non-dim]": sc_z,
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
        if self._solver_ == "Scikits":
            solver = pybamm.ScikitsDaeSolver(method="ida", atol=self.atoler, rtol=self.rtoler)
        elif self._solver_ == "Casadi":
            # model.convert_to_format = 'casadi'
            solver = pybamm.CasadiSolver(mode='fast', rtol=self.rtoler, atol=self.atoler, root_method="casadi")
    
        # Store discretised model and solver
        self._model = model
        self._param = param
        self._solver = solver

        #Store processed symbols/dimensionalization factors
        self._E_0 = param.process_symbol(self.E_0).evaluate()
        self._I_0 = param.process_symbol(self.I_0).evaluate()
        self._T_0 = param.process_symbol(self.T_0).evaluate()
        self._CS_d = param.process_symbol(self.CS_d).evaluate()
        # self._CP_d = param.process_symbol(self.CP_d).evaluate()
        # self._deltaT_nd = param.process_symbol(self.deltaT_nd).evaluate()

        # store time scale related things
        self._Tmax_nd = param.process_symbol(self.Tmax_nd).evaluate()
        self._m = self.t_steps
        self._x = self.x_steps
        
        print("Catalytic Model 08 initialized successfully.")

    def simulate(self, parameters):
        #####DEBUGGING#####
        # pybamm.set_logging_level("DEBUG")

        #7 May 23: method to pull times from init
        self.model("forward")
        times_nd = np.linspace(0, self._Tmax_nd, int(self._m)//2)
        times = np.linspace(0, self._Tmax_nd*2, int(self._m))
        print("Number of timesteps: " + str(self._m//2))
        print("Number of spacesteps: " + str(self._x))
        try:
            solution = self._solver.solve(self._model, times_nd, inputs=parameters)
            c_S = solution["S(soln) [non-dim]"](times_nd)
            c_Y = solution["Y(soln) [non-dim]"](times_nd)
            c_P = solution["P(ads) [non-dim]"](times_nd)
            c_X = solution["X(ads) [non-dim]"](times_nd)
            c_Z = solution["Z(ads) [non-dim]"](times_nd)
            cS = solution["S(soln) at electrode [non-dim]"](times_nd)
            E = solution["Applied Voltage [non-dim]"](times_nd)
            Ee = solution["Effective Voltage [non-dim]"](times_nd)
            current = solution["Current [non-dim]"](times_nd)
     
            self.model("backward", c_S[1:-1, -1], c_Y[1:-1, -1], c_P[-1], c_X[-1], c_Z[-1], Ee[-1])
            solution = self._solver.solve(self._model, times_nd, inputs=parameters)
            c_S = np.concatenate((c_S, solution["S(soln) [non-dim]"](times_nd)))
            c_P = np.concatenate((c_P, solution["P(ads) [non-dim]"](times_nd)))
            c_Z = np.concatenate((c_Z, solution["Z(ads) [non-dim]"](times_nd)))
            cS = np.concatenate((cS, solution["S(soln) at electrode [non-dim]"](times_nd)))
            E = np.concatenate((E, solution["Applied Voltage [non-dim]"](times_nd)))
            Ee = np.concatenate((Ee, solution["Effective Voltage [non-dim]"](times_nd)))
            current = np.concatenate((current, solution["Current [non-dim]"](times_nd)))
        
        except pybamm.SolverError as e:
            print(e)
            solution = np.zeros_like(times_nd)
            
        return (current, E, cS, c_P, times, times_nd)

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
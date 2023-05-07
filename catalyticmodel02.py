## This file defines the CatalyticModel class
#original author: Martin...?
#Model that has two normalizations, one wrt S

import pybamm
import numpy as np

class CatalyticModel:
    def __init__(self,const_parameters,seioptions):

        #Const_parameters are the parameters to be passed into this model via main.py
        param = pybamm.ParameterValues(const_parameters)

        #LE 31 Mar 23: idk if this is useful now
        #self.Npara = Npara

        #Options relating to the models, like SEI
        self.options= self.calloptions()

        # Create dimensional fixed parameters
        # the "_d" indicates the dimensional form, while the lack of one/prescence of "_nd" means it's nondimensional
        #Adding separate D parameters for S and P
        DS_d = pybamm.Parameter("Diffusion Coefficient of S [cm2 s-1]")
        DP_d = pybamm.Parameter("Diffusion Coefficient of P [cm2 s-1]")
        F = pybamm.Parameter("Faraday Constant [C mol-1]")
        R = pybamm.Parameter("Gas constant [J K-1 mol-1]")
        a = pybamm.Parameter("Electrode Area [cm2]")
        T = pybamm.Parameter("Temperature [K]")
        CS_d = pybamm.Parameter("Far-field concentration of S(soln) [mol cm-3]")
        CP_d = pybamm.Parameter("Far-field concentration of P(soln) [mol cm-3]")
        
        E_start_d = pybamm.Parameter("Voltage start [V]")
        E_reverse_d = pybamm.Parameter("Voltage reverse [V]")
        v = pybamm.Parameter("Scan Rate [V s-1]")
        Gamma = pybamm.Parameter("Electrode Coverage [mol cm-2]")

        # Create dimensional input parameters
        E0_d = pybamm.InputParameter("Reversible Potential [V]")
        k0_d = pybamm.InputParameter("Redox Rate [s-1]")
        kcat_forward_d = pybamm.InputParameter("Catalytic Rate For [cm3 mol-l s-1]")
        kcat_backward_d = pybamm.InputParameter("Catalytic Rate Back [cm3 mol-l s-1]")
        alpha = pybamm.InputParameter("Symmetry factor [non-dim]")
        Cdl_d = pybamm.InputParameter("Capacitance [F]")
        Ru_d = pybamm.InputParameter("Uncompensated Resistance [Ohm]")

        # Create scaling factors for non-dimensionalisation
        E_0 = (R * T)/ F #units are V
        T_0 = DS_d / a #units are seconds; CHANGED TO BE D/A
        Ctot = CS_d + CP_d #units are mol cm-3
        
        #get diffusion in here
        I_0 = (F * a * Gamma / T_0)  #units are A; this is (F^2) a Gammav / RT

        #creating time scale and non-dimensionalizing
        Tmax_d = abs(E_start_d - E_reverse_d)/v * 2
        Tmax_nd = Tmax_d / T_0
        
        # Non-dimensionalise parameters
        E0 = E0_d / E_0 #no units

        E_start = E_start_d / E_0
        E_reverse = E_reverse_d / E_0
        t_reverse = E_start - E_reverse

        k0 = k0_d * T_0 #no units
        kcat_for = kcat_forward_d * Ctot #no units
        kcat_back = kcat_backward_d * T_0 * Ctot #no units

        #Diffusion coefficients
        d_S = DS_d/DS_d #no units
        d_P = DP_d/DS_d #no units
        
        Cdl = Cdl_d * E_0 / (I_0 * T_0) #no units
        Ru = Ru_d * I_0 / E_0 #no units


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
        sc_Ox = pybamm.Variable("O(surf) [non-dim]")
        sc_Red = pybamm.Variable("R(surf) .[non-dim]")
        c_s = pybamm.Variable("S(soln) [non-dim]", domain="solution")
        c_p = pybamm.Variable("P(soln) [non-dim]", domain="solution")
        i = pybamm.Variable("Current [non-dim]")

        # Effective potential
        Eeff = Eapp - i * Ru #no units

        #"left" indicates environment directly on electrode surface; x = 0
        #"right" indicates environment between diffusion layer and bulk solution; x = xmax 

        # Faradaic current (Butler Volmer)
        i_f = k0 * ((sc_Red) * pybamm.exp((1-alpha) * (Eeff - E0)) #contribution of Red
                    - sc_Ox * pybamm.exp(-alpha * (Eeff - E0)) #contribution of Ox
                    )         

        # defining boundary values for S and P
        c_at_electrode_s = pybamm.BoundaryValue(c_s, "left")
        c_at_electrode_p = pybamm.BoundaryValue(c_p, "left")
        
        #catalytic rate contribution (this was previoulsly written as catalytic current)
        cat_con = kcat_for * c_at_electrode_s * (sc_Red) - kcat_back * c_at_electrode_p * (sc_Ox)

        # PDEs - left hand side is assumed to be time derivative of the PDE
        model.rhs = {
            sc_Ox: i_f + cat_con, #i_f is the echem contribution, cat_con is chemical contribution
            sc_Red: -i_f - cat_con,
            i: i_f - i, # capacitive current + Cdl * Eapp.diff(pybamm.t) 
            c_s: pybamm.div(pybamm.grad(c_s)),
            c_p: pybamm.div(pybamm.grad(c_p)),
        }

        # Setting boundary and initial conditions
        model.boundary_conditions = {
            c_s: {
                "right": (pybamm.Scalar(1), "Dirichlet"),
                "left": (-cat_con, "Neumann"),                    
            },

            c_p: {
                "right": (pybamm.Scalar(0), "Dirichlet"),   #0 makes sense - we'll always be starting with no product
                "left": (cat_con, "Neumann"),                 
            } 
        }

        model.initial_conditions = {
            sc_Ox: pybamm.Scalar(1),
            sc_Red: pybamm.Scalar(0),
            i: Cdl * Eapp.diff(pybamm.t), #again, capacitive current (if it's 0, starting i is 0)
            c_s: (CS_d/Ctot),
            c_p: (CP_d/Ctot),
        }

        # set spatial variables and domain geometry
        x = pybamm.SpatialVariable('x', domain="solution")
        x_max = 6*pybamm.sqrt(Tmax_nd)
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
            "O(surf) [non-dim]": sc_Ox,
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
        self._Ctot = param.process_symbol(Ctot).evaluate()

        # store time scale related things
        self._Tmax_d = param.process_symbol(Tmax_d).evaluate()
        
        print("Catalytic Model 02 initialized successfully.")

    def simulate(self, parameters):
        #####DEBUGGING#####
        #pybamm.set_logging_level("DEBUG")

        #7 May 23: method to pull times from init
        Times_d = np.linspace(0, self._Tmax_d, 2**12)
        times_nd = Times_d / self._T_0
        try:
            solution = self._solver.solve(self._model, times_nd, inputs=parameters)
        except pybamm.SolverError as e:
            print(e)
            solution = np.zeros_like(times_nd)
        return (
            solution["Current [non-dim]"](times_nd),
            solution["Applied Voltage [non-dim]"](times_nd),
            solution["O(surf) [non-dim]"](times_nd),
            solution["S(soln) at electrode [non-dim]"](times_nd),
            solution["P(soln) at electrode [non-dim]"](times_nd),
            times_nd
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
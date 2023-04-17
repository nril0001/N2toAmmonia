## This file defines the CatalyticModel class
#original author: Martin...?

import pybamm
import numpy as np

#This model builds off of the BaseBatteryModel in pybamm
class CatalyticModel:
    #initialization
    def __init__(self,const_parameters,seioptions):

        #Const_parameters are the parameters to be passed into this model via main.py
        param = pybamm.ParameterValues(const_parameters)

        #LE 31 Mar 23: idk if this is useful now
        #self.Npara = Npara

        #Options relating to the models, like SEI
        self.options= self.calloptions()

        # Create dimensional fixed parameters
        ##LE 20 Mar 23: D is now referenced in sigma (non dim scan rate)
        # the "_d" indicates the dimensional form, while the lack of one means it's nondimensional
        D = pybamm.Parameter("Diffusion Coefficient [cm2 s-1]")
        F = pybamm.Parameter("Faraday Constant [C mol-1]")
        R = pybamm.Parameter("Gas constant [J K-1 mol-1]")
        a = pybamm.Parameter("Electrode Area [cm2]")
        T = pybamm.Parameter("Temperature [K]")
        CS_d = pybamm.Parameter("Far-field concentration of S(soln) [mol cm-3]")
        CP_d = pybamm.Parameter("Far-field concentration of P(soln) [mol cm-3]")
        omega_d = pybamm.Parameter("Voltage frequency [rad s-1]")
        E_start_d = pybamm.Parameter("Voltage start [V]")
        E_reverse_d = pybamm.Parameter("Voltage reverse [V]")
        deltaE_d = pybamm.Parameter("Voltage amplitude [V]")
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
        E_0 = R * T / F #units are V
        T_0 = E_0 / v #units are seconds; this is RT/Fv
        I_0 = F * a * Gamma / T_0   #units are A; this is (F^2)aGammav / RT
        Ctot = CS_d + CP_d #units are mol cm-3
        L_0 = pybamm.sqrt(D * T_0) #units are cm
        #TODO:What is L_0 is meant to be a scalar for?

        # Non-dimensionalise parameters
        E0 = E0_d / E_0 #no units
    #old way of defining Es
        E_start = E_start_d / E_0
        E_reverse = E_reverse_d / E_0
        t_reverse = E_start - E_reverse
    #rewritten ways (with sigma)
        # E_start =  (E_start_d / E_0) - E0 #no units
        # E_reverse = (E_reverse_d / E_0) - E0 #no units
        # sigma = v * a / (E_0 * D) #no units
        # t_reverse = (E_start - E_reverse) / sigma #no units

        k0 = k0_d * T_0 #no units: D isn't in this because it's already on the surface
        kcat_for = kcat_forward_d * T_0 * Ctot #no units
        kcat_back = kcat_backward_d * T_0 * Ctot #no units
        
        Cdl = Cdl_d * E_0 / (I_0 * T_0) #21 Mar 21 - removed area, no units now
        Ru = Ru_d * I_0 / E_0 #no units


        #irrelevant for now
        # omega = 2 * np.pi * omega_d * T_0 #AC related
        # deltaE = deltaE_d / E_0 #no units, AC related

        # Input voltage protocol
        #pasted in from previous model...this works!
        #TODO: uses pybamm.t - understand how this variable functions better
        Edc_forward = -pybamm.t
        Edc_backwards = pybamm.t - 2 * t_reverse
        Eapp = E_start + \
            (pybamm.t <= t_reverse) * Edc_forward + \
            (pybamm.t > t_reverse) * Edc_backwards
            #deltaE * pybamm.sin(omega * pybamm.t)

        
    #2 Apr 23: Eapp rewritten to contain scanrate and nondim time
        # t_nd = pybamm.t / T_0
        # Edc_forward = E_start - (sigma * t_nd)
        # Edc_backwards = E_reverse + sigma * (t_nd - t_reverse)
        # Eapp = (t_nd <= t_reverse) * Edc_forward + \
        #     (t_nd > t_reverse) * Edc_backwards + \
        #     deltaE * pybamm.sin(omega * t_nd)                        

        # create PyBaMM model object
        model = pybamm.BaseBatteryModel(options=seioptions)
        #LE 21 Feb: switching from basemodel to basebatterymodel so that the SEI options actually matter

        ###Reaction of interest###
        # Ox* + e- <-> Red*
        # Red* + S -> Ox* + P

        ##ASSUMPTIONS FOR OX AND RED:
        # Both are completely confined on surface, so:
        #   No diffusion elements
        #   No spatial boundary conditions (only activity is at "left"/electrode surface)
        ###########################

        # Create state variables for model
        c_Ox = pybamm.Variable("O(surf) [non-dim]")
        c_s = pybamm.Variable("S(soln) [non-dim]", domain="solution")
        c_p = pybamm.Variable("P(soln) [non-dim]", domain="solution")
        i = pybamm.Variable("Current [non-dim]")

        #TODO: Why is this in i_f?
        # Effective potential
        Eeff = Eapp - i * Ru #no units

       
        #FIXME: check the boundary value for c_Ox is correct
        #"left" indicates environment directly on electrode surface; x = 0
        #"right" indicates environment between diffusion layer and bulk solution; x = xmax 

        # defining boundary values
        #c_Ox_boundary = pybamm.BoundaryValue(c_Ox, "left")
        c_at_electrode_s = pybamm.BoundaryValue(c_s, "left")
        c_at_electrode_p = pybamm.BoundaryValue(c_p, "left")

         # Faradaic current (Butler Volmer)
        #TODO: seems right, but this doesn't have indirect contributions from S or P
        i_f = k0 * ((1 - c_Ox) * pybamm.exp((1-alpha) * (Eeff - E0)) #contribution of Red
                    - c_Ox * pybamm.exp(-alpha * (Eeff - E0)) #contribution of Ox
                    )
        
        #catalytic rate contribution (this was previoulsly written as catalytic current)
        cat_con = kcat_for * c_at_electrode_s * (1-c_Ox) - kcat_back * c_at_electrode_p * (c_Ox)

        # PDEs - left hand side is assumed to be time derivative of the PDE
        model.rhs = {
            c_Ox: i_f + cat_con, #i_f is the echem contribution, cat_con is chemical contribution
            i: (i_f - cat_con - i)/(Cdl * Ru), # current divided by non-dim "s" obtained from Cdl and Ru term
            c_s: pybamm.div(pybamm.grad(c_s)), #TODO: understand this
            c_p: pybamm.div(pybamm.grad(c_p))
        }

        # algebraic equations (none)
        # model.algebraic = {
        # }

        # Setting boundary and initial conditions
        model.boundary_conditions = {
            c_Ox: {
                "right": (pybamm.Scalar(0), "Dirichlet"), #0 - no Ox should be present at xmax
                "left": (i_f, "Neumann"),
            },
            #LE 29 Mar 2023: copied c_s condition for c_Ox
            c_s: {
                "right": (pybamm.Scalar(1), "Dirichlet"),
                "left": (cat_con, "Neumann"),                    
            },

            c_p: {
                "right": (pybamm.Scalar(0), "Dirichlet"),   #0 makes sense - we'll always be starting with no product
                "left": (-cat_con, "Neumann"),                 
            } 
        }

        model.initial_conditions = {
            c_Ox: pybamm.Scalar(1),
            i: Cdl * Eapp.diff(pybamm.t), #having Cdl here is fine (if it's 0, starting i is 0)
            c_s: (CS_d/Ctot),
            c_p: (CP_d/Ctot),
        }

        # set spatial variables and domain geometry
        #LE 02 Apr 23: this can only handle one domain
        x = pybamm.SpatialVariable('x', domain="solution")
        model.geometry = pybamm.Geometry({
            "solution": {
                    x: {
                        "min": pybamm.Scalar(0),
                        "max": pybamm.Scalar(10)  #Add this in eventually: 6*sqrt(DTmax)
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
            "O(surf) [non-dim]": c_Ox,
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
        # self._omega_d = param.process_symbol(omega_d).evaluate()
        self._E_0 = param.process_symbol(E_0).evaluate()
        self._I_0 = param.process_symbol(I_0).evaluate()
        self._T_0 = param.process_symbol(T_0).evaluate()
        self._L_0 = param.process_symbol(L_0).evaluate()
        self._CS_d = param.process_symbol(CS_d).evaluate()
        self._CP_d = param.process_symbol(CP_d).evaluate()
        self._Ctot = param.process_symbol(Ctot).evaluate()

    #Actual simulation function
    def simulate(self, parameters, times):
        #####DEBUGGING#####
        #pybamm.set_logging_level("DEBUG")
        try:
            solution = self._solver.solve(self._model, times, inputs=parameters)
        except pybamm.SolverError as e:
            print(e)
            solution = np.zeros_like(times)
        return (
            solution["Current [non-dim]"](times),
            solution["Applied Voltage [non-dim]"](times),
            solution["O(surf) [non-dim]"](times),
            solution["S(soln) at electrode [non-dim]"](times),
            solution["P(soln) at electrode [non-dim]"](times),
            #parameters
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
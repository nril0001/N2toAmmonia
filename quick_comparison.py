## Main function
import catalyticmodel01 as cm1
import numpy as np
import time
import matplotlib.pylab as plt


def main(const_parameters,times):
     
    #FOR DIGIELCHM comparison
    files = [['digielchcomp/CV_kf_1_kb_1.txt',1, 1], ['digielchcomp/CV_kf_1e4_kb_4.1679.txt', 1e4, 4.1679],
             ['digielchcomp/CV_kf_1e5_kb_1.2578e-12.txt', 1e5, 1.2578e-12], ['digielchcomp/CV_kf_1e6_kb_4.4609e-20.txt', 1e6,4.4609e-20],
             ['digielchcomp/CV_kf_1e7_kb_1.5821e-27.txt', 1e7, 1.5821e-27], ['digielchcomp/CV_kf_100_kb_2.0416.txt', 100, 2.0416]]
    for i in files:
        
        #conditions that will change often over the course of testing
        input_parameters = {
            "Reversible Potential [V]": 0.0,
            "Redox Rate [s-1]": 10000,
            "Catalytic Rate For [cm3 mol-l s-1]": i[1],
            "Catalytic Rate Back [cm3 mol-l s-1]": i[2],
            "Symmetry factor [non-dim]": 0.5,
            #28 Mar 2023: not fully implemented
            "Capacitance [F]": 1e-8,
            "Uncompensated Resistance [Ohm]": 1.0
        }
        print(i)
        
        voltage = []
        curr = []
        row = []
         
        f = open(i[0],'r')
        for row in f:
            row = row.split("\t")
            voltage.append(float(row[0]))
            curr.append(float(row[1]))
                
        #list of options to pass into the model
        seioptions = ()
        #setting main model to reference CatalyticModel class
        cmodel = cm1.CatalyticModel(const_parameters,seioptions)
        #setting solved answers to ones usable here
        current, E_nd, O_nd, S_nd, P_nd = cmodel.simulate(input_parameters,times)
        ##redimensionalizing here for now. Messy to do in main, move later
        I_d = current * cmodel._I_0
        E_d = E_nd * cmodel._E_0
        title = i[0]
          
        
        plt.cla()
        plt.plot(voltage, curr, color = 'g', label = 'Digielch')
        plt.plot(E_d, I_d, color = 'b', label = 'Pybamm')
        plt.title(title)
        plt.legend()
        plt.xlabel("Eapp [V]")
        plt.ylabel("current [A]")
        plt.savefig(i[0]+".png")
    
    


if __name__ =='__main__':
    #Timer to measure performance
    ti = time.time()

    #constants that can vary, but generally won't change expt to expt
    const_parameters = {
        "Faraday Constant [C mol-1]": 96485.3328959,
        "Gas constant [J K-1 mol-1]": 8.314459848,
        "Far-field concentration of S(soln) [mol cm-3]": 1e-3,
        "Far-field concentration of P(soln) [mol cm-3]": 0e-6,
        "Diffusion Coefficient [cm2 s-1]": 1e-5,
        "Electrode Area [cm2]": 0.01,
        "Temperature [K]": 298,
        "Voltage frequency [rad s-1]": 9.0152,
        "Voltage start [V]": 0.5,
        "Voltage reverse [V]": -0.5,
        "Voltage amplitude [V]": 0.0,
        "Scan Rate [V s-1]": 0.05,
        "Electrode Coverage [mol cm-2]": 1e-9,
    }
    
    #BANDAID, FIXME#
    # calculate time scale to pass into model
    Tmax = abs(0.5 + 0.5)/0.05 * 2
    Tdim = np.linspace(0, Tmax, 2**12)
    TnonDim = (96485.3328959 * 0.05 / (8.314459848*298)) * Tdim

    

    main(const_parameters,TnonDim)
    

    #literally just to test that that main is working properly (delete later)
    print("complete in time: " + str((time.time()-ti)/60) + " minutes")

    #TODO: add multiple fxns for running models at the same time
    #TODO: improve plotting for pybamm based models: investigate using pybamm's plotting features??
    #TODO: LOW PRIORITY make compatible with interactive plottersz
## Main function
import attempt_at_adjusting as m
import numpy as np
import time
import matplotlib.pylab as plt


def main(const_parameters,times):
     
    #FOR DIGIELCHM comparison
    files = [['digielchcomp/current k0 1e-3.txt',1], ['digielchcomp/current k0 2e-3.txt', 2],
             ['digielchcomp/current k0 5e-3.txt', 5], ['digielchcomp/current k0_0.1.txt', 100],
             ['digielchcomp/current k0_0.5.txt', 500], ['digielchcomp/current k0_1.txt', 1000],
             ['digielchcomp/current k0_1e-2.txt', 10], ['digielchcomp/current k0_10.txt', 10000],
             ['digielchcomp/current k0_100.txt', 100000], ['digielchcomp/current k0_1000.txt', 1000000]]
    for i in files:
        
        #conditions that will change often over the course of testing
        input_parameters = {
            "Reversible Potential [V]": 0.5,
            "Redox Rate [s-1]": i[1],
            "Symmetry factor [non-dim]": 0.5,
            #28 Mar 2023: not fully implemented
            "Capacitance [F]": 1e-8,
            "Uncompensated Resistance [Ohm]": 1.0,
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
        #seioptions = ()
        #setting main model to reference CatalyticModel class
        cmodel = m.SingleReactionSolution(500, 6, const_parameters)
        #setting solved answers to ones usable here
        data = cmodel.simulate(input_parameters,times)
        ##redimensionalizing here for now. Messy to do in main, move later
        I_d = data['current']
        E_d = data['times'] * cmodel._E_0
        title = i[0]
          
        
        plt.cla()
        #plt.plot(voltage, curr, color = 'g', label = 'Digielch')
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
        "Far-field concentration of A [mol cm-3]": 1e-4,
        "Diffusion Coefficient [cm2 s-1]": 1e-5,
        "Electrode Area [cm2]": 1,
        "Temperature [K]": 298.2,
        "Voltage frequency [rad s-1]": 9.0152,
        "Voltage start [V]": 0.5,
        "Voltage reverse [V]": -0.5,
        "Voltage amplitude [V]": 0.0,
        "Scan Rate [V s-1]": 0.05,
        "Electrode Coverage [mol cm-2]": 1e-12,
    }
    
    #BANDAID, FIXME#
    # calculate time scale to pass into model
    n = 50
    t_eval = np.linspace(0, 50, n)
    
    Tmax = abs(0.5 + 0.5)/0.05 * 2
    Tdim = np.linspace(0, Tmax, 2**12)
    TnonDim = (96485.3328959 * 0.05 / (8.314459848*298)) * Tdim

    main(const_parameters,t_eval)
    

    #literally just to test that that main is working properly (delete later)
    print("complete in time: " + str((time.time()-ti)/60) + " minutes")

    #TODO: add multiple fxns for running models at the same time
    #TODO: improve plotting for pybamm based models: investigate using pybamm's plotting features??
    #TODO: LOW PRIORITY make compatible with interactive plottersz
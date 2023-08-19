## Main function
import catalyticmodel05_3_c as cm
import numpy as np
import time
import matplotlib.pylab as plt
import re
import os
from datetime import date


def main():
    #list of options to pass into the model
    seioptions = ()

    #retrieve date
    today = date.today()
    folder = today.strftime("%Y%m%d")

    #folder pathway - needs to be set to read all txt files that want analysed
    #use same naming conventions for text files as files in GammaTemp Variation
    pathway1 = "DigiElech/2023-08-15 small k0 no RC/CV/"
    pathway2 = "DigiElech/2023-08-15 small k0 no RC/SC/"
    output = "output/"+folder+"/"

    if os.path.isdir(output) == 0:
        os.mkdir(output, 0o666)
        os.mkdir(output+"SC/", 0o666)
        os.mkdir(output+"CV/", 0o666)

    #list to store files
    files_cv = []
    files_sc = []

    #Iterate directory of CVs
    for path in os.listdir(pathway1):
        #check if current path is a file
        if os.path.isfile(os.path.join(pathway1, path)):
            files_cv.append([pathway1+path])
            
    # Iterate directory of SCs
    for path in os.listdir(pathway2):
        #check if current path is a file
        if os.path.isfile(os.path.join(pathway2, path)):
            files_sc.append([pathway2+path])

    #retrieve variables from file pathway
    for t in files_cv:
        variables = re.findall(
            "[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", t[0])
        for l in range(len(variables)):
            variables[l] = float(variables[l])
            t.append(variables[l])
            
    x = [100]
    area = 1
    radius = np.sqrt(area/np.pi)
    Os = []
    Rs = []
    Es = []
    Is = []
    for o in x:
        ti = time.time()
        atol = 1e-12
        rtol = 1e-8
        t_steps = 2**14
        x_steps = o
        
        for i in files_cv:
            #constants that can vary, but generally won't change expt to expt
            print(i)
            const_parameters = {
                "Faraday Constant [C mol-1]": 96485.3328959,
                "Gas constant [J K-1 mol-1]": 8.314459848,
                "Far-field concentration of S(soln) [mol cm-3]": 1e-6,
                "Far-field concentration of P(soln) [mol cm-3]": 0,
                "Diffusion Coefficient of S [cm2 s-1]": i[9],
                "Diffusion Coefficient of P [cm2 s-1]": i[9],
                "Electrode Area [cm2]": 1,
                "Electrode Radius [cm]": radius,
                "Temperature [K]": 298,
                "Voltage start [V]": 0.5,
                "Voltage reverse [V]": -0.5,
                "Voltage amplitude [V]": 0.0,
                "Scan Rate [V s-1]": i[10],
                "Electrode Coverage [mol cm-2]": 1e-12,
            }
    
            #conditions that will change often over the course of testing     
            input_parameters = {
                "Reversible Potential [V]": 0.0,
                "Redox Rate [cm s-1]": i[6],
                "Catalytic Rate For [cm2 mol-l s-1]": 0,
                "Catalytic Rate Back [cm2 mol-l s-1]": 0,
                "Symmetry factor [non-dim]": 0.5,
                "Uncompensated Resistance [Ohm]": 1,
                "Capacitance [F]": 1,
            }
            
            # for unpacking DigiElech CVs
            voltage = []
            curr = []
            row = []
            count = 0
            f = open(i[0],'r')
            for row in f:
                count += 1
                if count < 3:
                    continue
                else:
                    row = row.split("\t")
                    voltage.append(float(row[0]))
                    curr.append(float(row[1]))
    
            #flipping current values from DigiElech
            curr = np.array(curr) 
    
            # for unpacking DigiElech SCs
            potential = []
            surfcon = []
            row = []
            count = 0
            f = open(files_sc[files_cv.index(i)][0],'r')
            for row in f:
                row = row.split("\t")
                try:
                    potential.append(float(row[0]))
                    surfcon.append(float(row[1]))
                except ValueError:
                    continue
                    
            # v = voltage
            # r = surfcon
    
            #setting main model to reference CatalyticModel class
            cmodel = cm.CatalyticModel(const_parameters,seioptions, atol, rtol, t_steps, x_steps)
            title = i[0].split("/")
            heading = files_sc[files_cv.index(i)][0].split("/")
            print(title[3])
            print(heading[3])
            current, E_nd, O_nd, R_nd, T_nd = cmodel.simulate(input_parameters)
            
            ##redimensionalizing here for now. Messy to do in main, move later
            I_d = current / cmodel._I_0
            E_d = E_nd / cmodel._E_0
            # Es.append(E_d)
            # Is.append(I_d)
            # Os.append(O_nd)
            # Rs.append(R_nd)
         
            #Plot current
            plt.cla()
            plt.plot(E_d[1:], (I_d[1:]), color = "red", label="PyBamm")
            plt.plot(voltage, np.array(-curr), color = 'blue', linestyle = 'dashdot', label = 'Digielch')
            plt.title("CV: K0_" + str(i[6]) + "_Ru_"+str(i[7]) +"_Cd_"+str(i[8]) + "_Ds_"+str(i[9])+"_Dp_"+str(i[9]) + "_SR_" + str(i[10]) )
            plt.xlabel("Potential [V]")
            plt.ylabel("Current [A]")
            plt.legend()
            plt.grid()
            plt.savefig(output+"CV/CV_cat05_dim_K0_" + str(i[6]) + "_Ru_"+str(i[7]) +"_Cd_"+str(i[8]) + "_Ds_"+str(i[9])+"_Dp_"+str(i[9]) + "_SR_" + str(i[10])+".png", dpi=600)
            np.savetxt(output+"CV/CV_cat05_dim_K0_" + str(i[6]) + "_Ru_"+str(i[7]) +"_Cd_"+str(i[8]) + "_Ds_"+str(i[9])+"_Dp_"+str(i[9]) + "_SR_" + str(i[10])+".dat", np.transpose(np.vstack((E_d, I_d))))
        
            # Plot concentration 
            plt.cla()
            plt.plot(E_d, (O_nd), color = "red", label="PyBamm - S")
            plt.plot(E_d, (R_nd), color = "orange", label="PyBamm - P")
            plt.plot(potential[:401], np.array(surfcon[:401])/surfcon[0], color = 'blue', linestyle = 'dashdot', label = 'Digielch - S')
            plt.plot(potential[401:], np.array(surfcon[401:])/surfcon[0], color = 'green', linestyle = 'dashdot', label = 'Digielch - P')
            plt.title("Norm SC: K0_" + str(i[6]) + "_Ru_"+str(i[7]) +"_Cd_"+str(i[8]) + "_Ds_"+str(i[9])+"_Dp_"+str(i[9]) + "_SR_" + str(i[10]) )
            plt.xlabel("Potential [V]")
            plt.ylabel("Surface Coverage [non-dim]")
            plt.legend()
            plt.grid()
            plt.savefig(output+"SC/SC_cat05_dim_K0_" + str(i[6]) + "_Ru_"+str(i[7]) +"_Cd_"+str(i[8]) + "_Ds_"+str(i[9])+"_Dp_"+str(i[9]) + "_SR_" + str(i[10])+".png", dpi=600)
            np.savetxt(output+"SC/SC_cat05_dim_K0_" + str(i[6]) + "_Ru_"+str(i[7]) +"_Cd_"+str(i[8]) + "_Ds_"+str(i[9])+"_Dp_"+str(i[9]) + "_SR_" + str(i[10])+".dat", np.transpose(np.vstack((E_d, O_nd, R_nd))))
            
            print("complete in time: " + str((time.time()-ti)/60) + " minutes")       
    
    #normalising the concentrations
    # plt.cla()
    # for u in range(len(Es)):    
    #     plt.plot(Es[u], (Os[u]), label=str(x[u]) + "- S")
    #     plt.plot(Es[u], (Os[u]), label=str(x[u]) + "- P")   
    # plt.plot(v[:398], np.array(r[:398])/r[0], color = 'blue', linestyle = 'dashdot', label = 'Digielch - S')
    # plt.plot(v[399:], np.array(r[399:])/r[0], color = 'green', linestyle = 'dashdot', label = 'Digielch - P')
    # plt.title("Normalised Concentrations of S at different X steps")
    # plt.xlabel("Potential [V]")
    # plt.ylabel("Surface Coverage [non-dim]")
    # plt.legend()
    # plt.grid()
    # plt.savefig(output+"Conc_S&P_cat05_dim_ko_"+str(i[5])+"_Ds_"+str(i[6])+"_Dp_"+str(i[7])+"._x_steps.png", dpi=600)    
    return

if __name__ =='__main__':
    #Timer to measure performance

    main()
    
    # print("complete in time: " + str((time.time()-ti)/60) + " minutes")

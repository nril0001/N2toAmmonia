## Main function
import catalyticmodel05 as cm
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
    pathway = "DigiElech/2023-06-06 Solution only/SC/"
    output = "output/"+folder+"/"

    if os.path.isdir(output) == 0:
        os.mkdir(output, 0o666)

    #list to store files
    files = []

    #Iterate directory
    for path in os.listdir(pathway):
        #check if current path is a file
        if os.path.isfile(os.path.join(pathway, path)):
            files.append([pathway+path])

    #retrieve variables from file pathway
    for t in files:
        variables = re.findall(
            "[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", t[0])
        for l in range(len(variables)):
            variables[l] = float(variables[l])
            t.append(variables[l])
    
    for i in files:
        #constants that can vary, but generally won't change expt to expt
        const_parameters = {
            "Faraday Constant [C mol-1]": 96485.3328959,
            "Gas constant [J K-1 mol-1]": 8.314459848,
            "Far-field concentration of S(soln) [mol cm-3]": 1e-6,
            "Far-field concentration of P(soln) [mol cm-3]": 0,
            "Diffusion Coefficient of S [cm2 s-1]": i[6],
            "Diffusion Coefficient of P [cm2 s-1]": i[7],
            "Electrode Area [cm2]": 1,
            "Temperature [K]": 298,
            "Voltage start [V]": 0.5,
            "Voltage reverse [V]": -0.5,
            "Voltage amplitude [V]": 0.0,
            "Scan Rate [V s-1]": 0.05,
            "Electrode Coverage [mol cm-2]": 1e-12,
        }

        #conditions that will change often over the course of testing     
        input_parameters = {
            "Reversible Potential [V]": 0.0,
            "Redox Rate [s-1]": i[5],
            "Catalytic Rate For [cm2 mol-l s-1]": 0,
            "Catalytic Rate Back [cm2 mol-l s-1]": 0,
            "Symmetry factor [non-dim]": 0.5,
            "Capacitance [F]": 0, #1e-8,
            "Uncompensated Resistance [Ohm]": 0
        }
        
        #for unpacking DigiElech CVs
        # voltage = []
        # curr = []
        # row = []
        # count = 0
        # f = open(i[0],'r')
        # for row in f:
        #     count += 1
        #     if count < 3:
        #         continue
        #     else:
        #         row = row.split("\t")
        #         voltage.append(float(row[0]))
        #         curr.append(float(row[1]))

        #flipping current values from DigiElech
        #curr = np.array(-curr)
        
        # potential = []
        # surf_cov = []
        # f = open(i[2],'r')
        # for row in f:
        #     row = row.split("\t")
        #     potential.append(float(row[0]))
        #     surf_cov.append(float(row[1])/100)

        #for unpacking DigiElech surface concentrations
        voltage = []
        surfcon = []
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
                surfcon.append(float(row[1]))

        #setting main model to reference CatalyticModel class
        cmodel = cm.CatalyticModel(const_parameters,seioptions)
        
        #setting solved answers to ones usable here
        #current, E_nd, O_nd, R_nd, S_nd, P_nd, cat_conc, i_f, k0, T_nd = cmodel.simulate(input_parameters)
        current, E_nd, O_nd, R_nd, T_nd = cmodel.simulate(input_parameters)
        
        ##redimensionalizing here for now. Messy to do in main, move later
        I_d = current * cmodel._I_0
        E_d = E_nd * cmodel._E_0

        # offset_E = []
        # maxx = E_d[np.where(I_d==max(I_d))[0][0]]
        # minn = E_d[np.where(I_d==min(I_d))[0][0]]
        # for l in range(0, 10001):
        #     offset_E.append(E_d[l] + maxx)
        # for l in range(10001, 20000):
        #     offset_E.append(E_d[l] + minn)
            
        temp = str(i[1])
        k = str(i[3])
        kcatf = str(i[4])
        kcatb = str(i[5])
        title = i[0].split("/")
        print(title[3])
        
    
        #QUICK PLOTS# 
        # plt.cla()
        # plt.plot(E_d, I_d, color = "b", label = "PyBamm")
        # #plt.plot(E_d, l*cmodel._I_0, color = "red", label = "test")
        # plt.plot(voltage, curr, color = 'g', label = 'Digielch')
        # plt.title(title[3])
        # plt.xlabel("Eapp [V]")
        # plt.ylabel("current [A]")
        # plt.grid()
        # plt.legend()
        # plt.savefig(output+"CV_cat05_dim_ko_"+str(i[5])+"_kf_"+str(i[6])+"_kb_"+str(i[7])+"_Cdl_"+str(i[9])+"_Ru_"+str(i[8])+"_Ds_"+str(i[10])+"_Dp_"+str(i[11])+".png", dpi=600)
        # np.savetxt(output+"CV_cat04_dim_k0_.dat", np.transpose(np.vstack((E_d, S_nd))))
        
        #normalising the concentrations
        maxval = np.max(np.append(O_nd,R_nd))
        plt.cla()
        plt.plot(E_d, (O_nd/maxval), color = "red", label="PyBamm - S")
        plt.plot(E_d, (R_nd/maxval), color = "orange", label="PyBamm - P")
        plt.plot(voltage[:398], np.array(surfcon[:398])/np.max(surfcon), color = 'blue', linestyle = 'dashdot', label = 'Digielch - S')
        plt.plot(voltage[399:], np.array(surfcon[399:])/np.max(surfcon), color = 'green', linestyle = 'dashdot', label = 'Digielch - P')
        #plt.plot(E_d, cat, color = "orange", label="PyBamm - Cat")
        plt.title("Normalised Concentrations of S and P")
        plt.xlabel("Potential [V]")
        plt.ylabel("Surface Coverage [non-dim]")
        plt.legend()
        plt.grid()
        #plt.savefig(output+"Conc_S&P_cat05_dim_ko_"+str(i[5])+"_kf_"+str(i[6])+"_kb_"+str(i[7])+"_Cdl_"+str(i[9])+"_Ru_"+str(i[8])+"_Ds_"+str(i[10])+"_Dp_"+str(i[11])+".png", dpi=600)
        plt.savefig(output+"Conc_S&P_cat05_dim_ko_"+str(i[5])+"_Ds_"+str(i[6])+"_Dp_"+str(i[7])+".png", dpi=600)
        
        # plt.cla()
        # plt.plot(offset_E[:10001], (I_d[:10001]), color = "red", label = "Pybamm Maxmium peak (Oxidation)")
        # plt.plot(offset_E[10001:], (I_d[10001:]), color = "orange", label = "Pybamm Minimum peak (Reduction)")
        # plt.plot(voltage, curr, color = 'g', label = 'Digielch')
        # plt.title("Offset Voltammogram k0: " + k00)
        # plt.xlabel("Eapp [V]")
        # plt.ylabel("current [A]")
        # plt.grid()
        # plt.legend()
        # plt.savefig(output+"CV_cat04_Offset_k0_"+k00+".png", dpi=600)
        
        # plt.cla()
        # plt.plot(offset_E[:10001], abs(I_d[:10001]), color = "red", label = "Maxmium peak (Oxidation)")
        # plt.plot(offset_E[10001:], abs(I_d[10001:]), color = "orange", label = "Minimum peak (Reduction)")
        # plt.title("k0: " + k00)
        # plt.xlabel("Eapp [V]")
        # plt.ylabel("current [A]")
        # plt.grid()
        # plt.legend()
        # plt.savefig(output+"CV_cat04_offset_k0_"+k00+".png", dpi=600)
        # np.savetxt(output+"CV_cat04_offset_k0_"+k00+".dat", np.transpose(np.vstack((E_d, I_d))))
    
        # plt.cla()
        # plt.plot(E_nd, current)
        # plt.xlabel("Eapp [non-dim]")
        # plt.ylabel("current [non-dim]")
        # plt.grid()
        # plt.savefig(output+"CV_cat04_nondim_cat01.png", dpi=600)
        # # np.savetxt(output+"CV_cat04_nondim_pybamm_kf_1.dat", np.transpose(np.vstack((E_nd, current))))
    
        # plt.cla()
        # plt.plot(T_nd, current)
        # plt.xlabel("time [non-dim]")
        # plt.ylabel("current [non-dim]")
        # plt.grid()
        # plt.savefig(output+"currentvstime_cat03.png", dpi=600)
        
        # plt.cla()
        # plt.plot(T_nd, E_nd)
        # plt.xlabel("time [non-dim]")
        # plt.ylabel("Eapp [non-dim]")
        # plt.grid()
        # plt.savefig(output+"Eappvstime_cat03.png", dpi=600)
        
        # plt.cla()
        # plt.plot(T_nd, cat_conc, label="Cat_conc")
        # plt.plot(T_nd, i_f, label="i_f")
        # plt.xlabel("time [non-dim]")
        # plt.ylabel("rates [non-dim]")
        # plt.legend()
        # plt.grid()
        # plt.savefig(output+"ratesvstime_cat03.png", dpi=600)
        # # np.savetxt(output+"ratesvstime_cat01.dat", np.transpose(np.vstack((T_nd, O_nd, R_nd))))
        
        # plt.cla()
        # #plt.plot(voltage, curr, color = 'g', label = 'Digielch')
        # plt.plot(E_d, I_d, color = 'b', label = 'Pybamm')
        # #plt.title(title)
        # plt.legend()
        # plt.xlabel("Eapp [V]")
        # plt.ylabel("current [A]")
        # plt.savefig(output+i[0]+".png")
        
    return

if __name__ =='__main__':
    #Timer to measure performance
    ti = time.time()

    main()
    
    #literally just to test that that main is working properly (delete later)
    print("complete in time: " + str((time.time()-ti)/60) + " minutes")

    #TODO: add multiple fxns for running models at the same time
    #TODO: improve plotting for pybamm based models: investigate using pybamm's plotting features??
    #TODO: LOW PRIORITY make compatible with interactive plottersz
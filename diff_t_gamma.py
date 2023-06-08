# -*- coding: utf-8 -*-
"""
Created on Wed May 24 10:55:34 2023

@author: Nathan
"""

## Main function
import catalyticmodel04 as cm
import numpy as np
import time
import matplotlib.pylab as plt
import os
import re
from datetime import date

def main():
    #list of options to pass into the model
    seioptions = ()
    
    #retrieve date
    today = date.today()
    folder = today.strftime("%Y%m%d")
    
    #folder pathway - needs to be set to read all txt files that want analysed
    #use same naming conventions for text files as files in GammaTemp Variation
    pathway = "DigiElech/surfaceonly/GammaTemp Variation/"
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
        variables = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", t[0])
        for l in range(len(variables)):
            variables[l] = float(variables[l])
            t.append(variables[l])  
    
    for i in files:
        const_parameters = {
            "Faraday Constant [C mol-1]": 96485.3328959,
            "Gas constant [J K-1 mol-1]": 8.314459848,
            "Far-field concentration of S(soln) [mol cm-3]": 1e-2,
            "Far-field concentration of P(soln) [mol cm-3]": 0e-6,
            "Diffusion Coefficient of S [cm2 s-1]": 1e-5,
            "Diffusion Coefficient of P [cm2 s-1]": 1e-5,
            "Electrode Area [cm2]": 1,
            "Temperature [K]": i[1],
            "Voltage start [V]": 0.5,
            "Voltage reverse [V]": -0.5,
            "Voltage amplitude [V]": 0.0,
            "Scan Rate [V s-1]": 0.05,
            "Electrode Coverage [mol cm-2]": i[2],
        }
        
        
        input_parameters = {
            "Reversible Potential [V]": 0.0,
            "Redox Rate [s-1]": 1000,
            "Catalytic Rate For [cm2 mol-l s-1]": 0,
            "Catalytic Rate Back [cm2 mol-l s-1]": 0,
            "Symmetry factor [non-dim]": 0.5,
            "Capacitance [F]": 0, #1e-8,
            "Uncompensated Resistance [Ohm]": 0.0
        }
        
        # read potential and current from digielch files
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
        curr = np.array(curr) * (-1)
        
        #digielch surface concentration profiles
        # potential = []
        # surf_cov = []
        # f = open(i[0],'r')
        # for row in f[2:]:
        #     row = row.split("\t")
        #     potential.append(float(row[0]))
        #     surf_cov.append(float(row[1])/100)
        
        
        #setting main model to reference CatalyticModel class
        cmodel = cm.CatalyticModel(const_parameters,seioptions)
        #setting solved answers to ones usable here
        current, E_nd, O_nd, R_nd, S_nd, P_nd, cat_conc, i_f, k0, T_nd = cmodel.simulate(input_parameters)
        ##redimensionalizing here for now. Messy to do in main, move later
        I_d = current * cmodel._I_0
        print(cmodel._I_0)
        #I_0 = (F * a * Gamma * v) / ((R * T)/F)  #units are A
        E_d = E_nd * cmodel._E_0

        offset_E = []
        maxx = E_d[np.where(I_d==max(I_d))[0][0]]
        minn = E_d[np.where(I_d==min(I_d))[0][0]]
        for l in range(0, 10001):
            offset_E.append(E_d[l] + maxx)
        for l in range(10001, 20000):
            offset_E.append(E_d[l] + minn)
            
        k00 = str(i[1])
        k01 = str(i[2])
    
        #QUICK PLOTS# 
        plt.cla()
        plt.plot(E_d, I_d, color = "b", label = "PyBamm")
        plt.plot(voltage, curr, color = 'g', label = 'Digielch')
        plt.title("Voltammogram k0: " + k00)
        plt.xlabel("Eapp [V]")
        plt.ylabel("current [A]")
        plt.grid()
        plt.legend()
        plt.savefig(output+"CV_dim_"+k00+"K_G_"+k01+".png", dpi=1200)
        #np.savetxt(output+"CV_dim_"+k00+"K_G_"+k01+".dat", np.transpose(np.vstack((E_d, I_d))))
        
        # plt.cla()
        # plt.plot(offset_E[:10001], (I_d[:10001]), color = "red", label = "Pybamm Maxmium peak (Oxidation)")
        # plt.plot(offset_E[10001:], (I_d[10001:]), color = "orange", label = "Pybamm Minimum peak (Reduction)")
        # plt.plot(voltage, curr, color = 'g', label = 'Digielch')
        # plt.title("Offset Voltammogram k0: " + k00)
        # plt.xlabel("Eapp [V]")
        # plt.ylabel("current [A]")
        # plt.grid()
        # plt.legend()
        # plt.savefig(output+"CV_offset_dim_"+k00+"K_G_"+k01+".png", dpi=600)
        
        # plt.cla()
        # plt.plot(offset_E[:10001], abs(I_d[:10001]), color = "red", label = "Maxmium peak (Oxidation)")
        # plt.plot(offset_E[10001:], abs(I_d[10001:]), color = "orange", label = "Minimum peak (Reduction)")
        # plt.title("k0: " + k00)
        # plt.xlabel("Eapp [V]")
        # plt.ylabel("current [A]")
        # plt.grid()
        # plt.legend()
        # plt.savefig(output+"CV_offset_flipped_dim_"+k00+"K_G_"+k01+".png", dpi=600)
        # np.savetxt(output+"CV_offset_flipped_dim_"+k00+"K_G_"+k01+".dat", np.transpose(np.vstack((E_d, I_d))))
        
        # plt.cla()
        # plt.plot(E_d, O_nd, color = "red", label="PyBamm - Ox")
        # plt.plot(E_d, R_nd, color = "orange", label="PyBamm - Red")
        # plt.plot(potential[:20000], surf_cov[:20000], color = "blue", label="DigiElch - Ox")
        # plt.plot(potential[20000:], surf_cov[20000:], color = "green", label="DigiElch - Red")
        # plt.title("Concentration Profile: " + k00)
        # plt.xlabel("Potential [V]")
        # plt.ylabel("Surface Coverage [non-dim]")
        # plt.legend()
        # plt.grid()
        # plt.savefig(output+"Surf_cov_dim_"+k00+"K_G_"+k01+".png", dpi=600)
    

        # plt.cla()
        # plt.plot(E_nd, current)
        # plt.xlabel("Eapp [non-dim]")
        # plt.ylabel("current [non-dim]")
        # plt.grid()
        # plt.savefig(output+"currentvsEapp_cat01.png", dpi=600)
        # # np.savetxt(output+"current_nondim_pybamm_kf_1.dat", np.transpose(np.vstack((E_nd, current))))
    
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
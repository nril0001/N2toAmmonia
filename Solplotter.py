# -*- coding: utf-8 -*-
"""
Created on Wed May  6 02:12:26 2020

@author: luke
"""
import numpy as np
import scipy as sci
import Sollplotter_mod as mod
import matplotlib.pyplot as plt

def AnalyticalModel(Kcat1):
# Parameters (Use volts)
    Estart = -0.5
    Eend = 0.5
    E0 = 0
    scanrate = 0.05
    Kcat1 = Kcat1 # (s^-1)
    dE = 0.0 #Sine wave ampliture
    freq = 9
    
    # constants
    F = 96485.3328959 #faradays constants
    Rg = 8.31447  #Gas constant
    T = 298.2    #temperature 
    
    # odd stuff
    res = 1000
    trun = 10
    
    Et = np.linspace(Estart,Eend,res)   # non-normalised scanrate
    
    # Nondimensionalising parameters
    dE = F/(Rg*T)*dE
    Eadjust = (F/(Rg*T))*(Et - E0)/2
    Kcat1 =((Rg*T)/(F*scanrate))*Kcat1
    Nondimang = ((Rg*T)/(F*scanrate))*(2*np.pi*freq)
    
    
    
    # calculate IDCNorm
    IDCnorm = np.ones(res)*(-Kcat1/2)
    
    for N in range(0,trun+1):
        # case for when g is just tanh
        if N == 0:
            g1 = np.tanh(Eadjust) 
        else:
            g1 = mod.Tanh_nthdev(2*N,Eadjust) 
        g2 = mod.Tanh_nthdev(2*N+1,Eadjust)
        """NEED TO DOUBLE CHECK IF ITS 2 TIMES N FACTORIAL OR 2N factorial"""
        holder = ((dE/4)**(2*N))*(g2- 2*Kcat1*g1)/((2*sci.special.factorial(N))**2)
        IDCnorm  -= holder
    
    #  # Plot odd harmonics
    # Ioddnorm = np.zeros(res)
    # for N in range(0, trun+1):
    #     g2 = mod.Tanh_nthdev(2*N+1, Eadjust)*((dE/4)**(2*N+1))
    #     holder = 0
    #     for i in range(0, N+1):
    #         scalar = ((-1) ^ i)/(sci.special.factorial(N-i)
    #                     * sci.special.factorial(N+i+1))
    #         AC = Kcat1*np.sin((2*i+1)*omega*TnonDim) + (2*i+1) * \
    #                             omega*np.cos((2*i+1)*omega*TnonDim)
    #         holder += scalar*AC
    #     Ioddnorm += g2*holder
    # # Plot even harmonics
    # Ievennorm = np.zeros(res)
    # for N in range(0, trun+1):
    #     g2 = mod.Tanh_nthdev(2*N+2, Eadjust)*((dE/4)**(2*N+2))
    #     holder = 0
    #     for i in range(0, N+1):
    #         scalar = ((-1) ^ i)/(sci.special.factorial(N-i)
    #                     * sci.special.factorial(N+i+2))
    #         AC = -Kcat1*np.cos((2*i+2)*omega*TnonDim) + (2*i+2) * \
    #                             omega*np.sin((2*i+2)*omega*TnonDim)
    #         holder += scalar*AC
    #     Ievennorm += g2*holder
    # # extract non dim harmonics - These are all in nondimensional currents
    # bandwidth = np.array([0, 4, 0, 4, 0, 4, 0, 4, 0])
    # Harmonicsodd = mod.harm_gen(Ioddnorm, Tdim, [freq], bandwidth, 1)
    
    # bandwidth = np.array([0, 0, 4, 0, 4, 0, 4, 0, 4])
    # Harmonicseven = mod.harm_gen(Ievennorm, Tdim, [freq], bandwidth, 1)
    # Harmonicseven = mod.harm_gen(Ievennorm,TnonDim,[omega],bandwidth, 1)
    return Et, IDCnorm
    #plt.plot(Et,IDCnorm) #+Ioddnorm+Ievennorm)
    #plt.savefig("output/solplot.png")
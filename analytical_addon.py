"""
This code is dependent on recieving a voltage protocol from pybamm or another source, so it's more like an add-on rather than a standalone model.
Based on analytical solution from J. Zhang, A.M. Bond / Journal of Electroanalytical Chemistry 600 (2007) 23-34

Original author: Luke Gundry
"""
import Sollplotter_mod as mod
import numpy as np
import scipy as sci

class AnalyticalModel:
    #pass in Eadjust and Kcat1 (dimensionless)
    def simulate(Eadjust,Kcat1):
        res = 2**12  # Needs to be high enough or fft doesn't go to high enough frequency
        dE = 0 #sine wave amplitude, 0 for DC
        # odd stuff
        trun = 8  # DON"T INCREASE PAST THIS VALUE ITS PRETTY UNSTABLE

        # calculate IDCNorm
        IDCnorm = np.ones(res)*(-Kcat1/2)

        for N in range(0, trun+1):
            # case for when g is just tanh
            if N == 0:
                g1 = np.tanh(Eadjust)
            else:
                g1 = mod.Tanh_nthdev(2*N, Eadjust)
            g2 = mod.Tanh_nthdev(2*N+1, Eadjust)
            holder = ((dE/4)**(2*N))*(g2 - 2*Kcat1*g1)/((2*sci.special.factorial(N))**2)
            IDCnorm -= holder
        #------------------------------------------------------------------
        ##LE 29 Jan 23: AC contributions    
        # # Plot odd harmonics
        # Ioddnorm = np.zeros(res)
        # for N in range(0, trun+1):
        #     g2 = mod.Tanh_nthdev(2*N+1, Eadjust)*((dE/4)**(2*N+1))
        #     holder = 0
        #     for i in range(0, N+1):
        #         scalar = ((-1) ^ i)/(sci.special.factorial(N-i)
        #                   * sci.special.factorial(N+i+1))
        #         AC = Kcat1*np.sin((2*i+1)*omega*TnonDim) + (2*i+1) * \
        #                           omega*np.cos((2*i+1)*omega*TnonDim)
        #         holder += scalar*AC
        #     Ioddnorm += g2*holder
        # # Plot even harmonics
        # Ievennorm = np.zeros(res)
        # for N in range(0, trun+1):
        #     g2 = mod.Tanh_nthdev(2*N+2, Eadjust)*((dE/4)**(2*N+2))
        #     holder = 0
        #     for i in range(0, N+1):
        #         scalar = ((-1) ^ i)/(sci.special.factorial(N-i)
        #                   * sci.special.factorial(N+i+2))
        #         AC = -Kcat1*np.cos((2*i+2)*omega*TnonDim) + (2*i+2) * \
        #                            omega*np.sin((2*i+2)*omega*TnonDim)
        #         holder += scalar*AC
        #     Ievennorm += g2*holder
        # # extract non dim harmonics - These are all in nondimensional currents
        # bandwidth = np.array([0, 4, 0, 4, 0, 4, 0, 4, 0])
        # Harmonicsodd = mod.harm_gen(Ioddnorm, Tdim, [freq], bandwidth, 1)

        # bandwidth = np.array([0, 0, 4, 0, 4, 0, 4, 0, 4])
        # Harmonicseven = mod.harm_gen(Ievennorm, Tdim, [freq], bandwidth, 1)
        # Harmonicseven = mod.harm_gen(Ievennorm,TnonDim,[omega],bandwidth, 1)
        ##LE 30 Jan 23: Harmonics stuff over
        #------------------------------------------------------------
        return (IDCnorm)
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 02:12:53 2020

@author: luke
"""
import scipy as sci
from scipy.special import binom
from numpy import exp,tanh,zeros,ones
import numpy as np
from scipy.fftpack import rfft, irfft, rfftfreq
from scipy.signal import hilbert as anal_hil_tran # NEED signal to make a smooth envolope
import matplotlib.pyplot as plt
# gives the nth order tanh derivative Via a known function
""""NOTE THAT DUE TO THE DIFFERENCE IN THESE EQUATIONS THERE IS DIFFERENCE 
BETWEEN TAN AND EXP SOLUTION OF ORDER OF 1*10^-11 %"""
def Tanh_nthdev(order,array):
    # does the scalar thingy
    top =2**(order+1)*exp(2*array)
    bottom = (1 + exp(2*array))**(order + 1)
    scalar = top/bottom
    
    # sum part
    x  = zeros(len(array))
    for k in range(0,order):
        eul = eulerian_num(order,k)
        x[:] += ((-1)**k)*eul*exp(2*k*array[:])
        
    return scalar*x     # returns scalar times sum

# gives the eulerian number
def eulerian_num(n,k):
    x = 0
    for j in range(0,k+2):
        c = binom(n+1,j)
        x += (-1)**j *(k+1-j)**n*c  
    return x

def ifft_windowing(AC_freq, freq, bandwidth,Nac):
    
    nharm = len(bandwidth) - 1
    g = freq[1] # gets gradent for translation to n space 
    f_window = np.zeros((Nac*nharm, 2))
    
    j = 0   #  AC counter
    while j != Nac:   # collects the AC input freqency
        f = AC_freq[j]
        
        i = 0   # harmonic counter
        while i != nharm:   # collects the harmonics for each freqency
            width = bandwidth[i + 1]/2
            f_window[(j*nharm + i),0] = (i + 1)*f - width
            f_window[(j*nharm + i),1] = (i + 1)*f + width
            i += 1
        j += 1
    
    
    # move frequency to n space
    n_window = (2/g)*f_window         #ROUND COULD BE THERE
    
    return n_window, f_window, g

def Harmonic_gen(fft_res, n_window, Nac, bandwidth, g): # cuts the harmonic fourier space data out {NEEDS TO EXCLUDE FUNDEMENTAL}
    
    nharm = len(bandwidth) - 1     # counts the number of harmonics
    Np = len(fft_res)    # gets the number of datapoints
    # harmstore = np.zeros(nharm*spaces[4] ,int(Mwin*(2/g)))
    hil_store = np.zeros((nharm*Nac + 1 ,Np))
    N_ifft = np.zeros((nharm*Nac + 1))
    
    
    #extracts the fundimental haronic
    N_fund = int(round(bandwidth[0]/g))
    x  = fft_res[0:N_fund]  # Cuts out windowed harmonics 
    y = np.zeros(Np)
    jj=0
    while jj != N_fund:  
        y[jj] = x[jj]
        jj += 1
    hil_store[0,:] = irfft(y)  # generates fundimental harmonics from cut window
       
    
    j = 0
    while j != Nac:
    # something to take fft and cut windows into the four space thing
        
        i = 0
        while i!= nharm:
            wl = int(round(n_window[ i,0]))
            wh = int(round(n_window[i,1]))
            x = fft_res[wl:wh]   #This need to be truncated in the future
            y = np.zeros(Np)
            jj=wl
            while jj != wh:
                y[jj] = x[jj-wl]
                jj += 1
            harmonic = irfft(y)   # generates harmonics
            
            hil_store[(j*nharm + i) + 1, :] = abs(anal_hil_tran(harmonic)) #uses HILBERT TRANSFORM to generate the envolope
            
            # using the abs fixed an issue with the complexes disapearing in the return
            
            i += 1
        j += 1
 
    return hil_store, N_ifft

def harm_gen(Curr,time,AC_freq,bandwidth, Nac):
    
    fft_res = rfft(Curr)
    # this can be done seperatly So lok into it
    freq = rfftfreq(len(fft_res),d = time[1])   # for ML could be done at the start to stop iterations ALSO DOUBLE VECTOR
    #plt.plot(freq[0:5000],np.log10(fft_res[0:5000]).real)
    n_window, f_window, g =  ifft_windowing(AC_freq, freq, bandwidth, Nac)     # Extracts the frequency and n windowing {CAN BE SEPERATE}
    
    hil_store, N_ifft = Harmonic_gen(fft_res, n_window,Nac, bandwidth, g)  # extract envolope data
    
    return hil_store

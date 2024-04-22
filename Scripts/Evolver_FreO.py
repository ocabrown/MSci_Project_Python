#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 16:46:31 2022

@author: ollie
"""

import L_Finder_FreO as LFFO
import Maximum_Accretion_Rate_Calculator as MARC
import numpy as np
import scipy.constants as spc
import matplotlib.pyplot as plt




def Evolver(mc,me,rp,pd,kd,Td,L_guess_ini,Z,acc_k,rc,acc,n_i,m_grid,t_tot,n_t,Sigma_in,Ms,Rp,r_in,r_out,alpha,num_sig):
    
    L0, vari, _, _, _, _ = LFFO.L_Finder(mc,me,rp,pd,kd,Td,L_guess_ini,Z,acc_k,rc,acc,n_i,m_grid, Plot=False)
    
    n_t += 1
    
    empty = np.zeros((n_t))
    empty_vari = np.array(vari) * 0
    xy = np.shape(empty_vari)
    
    R = rc                          # Core Radius [cm]
    G = spc.G * 1e3                 # Gravitational constant [cm^3/g/s^2]
    M_E = 5.972e27                  # Earth mass [g]
    L_Sun = 3.846e33                # Solar luminosity [erg/s]
    t_tot *= (365*24*60*60)         # Convert input total time in years into s
    dt = t_tot / (n_t-1)
    
    M, L, t = np.array([empty, empty, empty])
    varis = np.full((n_t,xy[0],xy[1]), empty_vari)
    M[0] = (mc + me) * M_E
    L[0], varis[0] = L0, vari
    tn = 0
    
    plt.figure()
    
    for n in range(n_t-1):
        
        Ln = L[n] * L_Sun
        M[n+1] = M[n] + (Ln*R*dt)/(G*M[n])
        
        Mdot_max, Sigma, R_sig = MARC.Max_Acc_Rate(Sigma_in,Ms,M[n],Td,Rp,r_in,r_out,alpha,num_sig)
        factor = 365*24*60*60/M_E
        print(M[n]/M_E,((M[n+1]-M[n])/dt)*factor,Mdot_max*factor)
        
        if n == 0:
            plt.plot(R_sig, Sigma, color="blue", label=f"{t[n]/(1e6*365*24*60*60)} Myr")
            plt.axvline(R_sig[int(num_sig/5)], color="blue", linestyle="--")
            plt.axvline(R_sig[int(2*num_sig/5)-1], color="blue", linestyle="--")
            plt.axvline(R_sig[int(3*num_sig/5)-1], color="blue", linestyle="--")
            plt.axvline(R_sig[int(4*num_sig/5)-1], color="blue", linestyle="--")
            
        if n == 3:
            plt.plot(R_sig, Sigma, color="orange", label=f"{t[n]/(1e6*365*24*60*60)} Myr")
            plt.axvline(R_sig[int(num_sig/5)], color="orange", linestyle="--")
            plt.axvline(R_sig[int(2*num_sig/5)-1], color="orange", linestyle="--")
            plt.axvline(R_sig[int(3*num_sig/5)-1], color="orange", linestyle="--")
            plt.axvline(R_sig[int(4*num_sig/5)-1], color="orange", linestyle="--")

        if n == 6:
            plt.plot(R_sig, Sigma, color="green", label=f"{t[n]/(1e6*365*24*60*60)} Myr")
            plt.axvline(R_sig[int(num_sig/5)], color="green", linestyle="--")
            plt.axvline(R_sig[int(2*num_sig/5)-1], color="green", linestyle="--")
            plt.axvline(R_sig[int(3*num_sig/5)-1], color="green", linestyle="--")
            plt.axvline(R_sig[int(4*num_sig/5)-1], color="green", linestyle="--")
            
        if n == 9:
            plt.plot(R_sig, Sigma, color="red", label=f"{t[n]/(1e6*365*24*60*60)} Myr")
            plt.axvline(R_sig[int(num_sig/5)], color="red", linestyle="--")
            plt.axvline(R_sig[int(2*num_sig/5)-1], color="red", linestyle="--")
            plt.axvline(R_sig[int(3*num_sig/5)-1], color="red", linestyle="--")
            plt.axvline(R_sig[int(4*num_sig/5)-1], color="red", linestyle="--")
        
        if ((M[n+1]-M[n])/dt) > Mdot_max:
            print("Max Reached!")
            M[n+1] = M[n] + (dt * Mdot_max)
        
        me += (M[n+1] - M[n]) / M_E
        
        m_grid *= M[n+1]/M[n]
        
        L[n+1], varis[n+1], _, _, _, _ = LFFO.L_Finder(mc,me,rp,pd,kd,Td,L[n],Z,acc_k,rc,acc,n_i,m_grid, Plot=False)
        
        tn += dt
        t[n+1] = tn
    
    plt.xlabel("Radius / AU")
    plt.ylabel("$\Sigma$ / $gcm^{-2}$")
    plt.legend()
    plt.show()
    
    t /= (365*24*60*60)             # Convert time in s into years for plotting
    
    return L, M, t, varis
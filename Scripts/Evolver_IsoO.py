#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 16:46:31 2022

@author: ollie
"""

import Mass_Grid_Solver as MGS
import L_Finder_IsoO as LFIO
import Maximum_Accretion_Rate_Calculator as MARC
import numpy as np
import scipy.constants as spc
import matplotlib.pyplot as plt





def EvolverIsoO(mc,me,rp,pd,kd,Td,L_guess_ini,acc_l,n_i,acc_m,p_ana_c_ini,m_grid,t_tot,n_t,Sigma_in,Ms,Rp,r_in,r_out,alpha,num_sig,L09=False):
    
    #m_grid, p_ana_c_ini = MGS.Mass_Grid_Solver_Iso_Opacity(mc,me,rp,pd,kd,Td,L_guess_ini,n_i,acc_m,p_ana_c_ini)
    L0, vari = LFIO.L_Finder(mc,me,rp,pd,kd,Td,L_guess_ini,acc_l,n_i,m_grid, Plot=False)[:2]
    
    n_t += 1
    
    empty = np.zeros((n_t))
    empty_vari = np.array(vari) * 0
    xy = np.shape(empty_vari)
    
    M_E = 5.972e27                      # Earth mass [g]
    pc = 3.2                            # Core density [g/cm^3]
    R = (3*mc*M_E/(4*np.pi*pc))**(1/3)  # Core Radius [cm]
    
    G = spc.G * 1e3                 # Gravitational constant [cm^3/g/s^2]
    M_E = 5.972e27                  # Earth mass [g]
    L_Sun = 3.846e33                # Solar luminosity [erg/s]
    t_tot *= (365*24*60*60)         # Convert input total time in years into s
    dt = t_tot / (n_t-1)
    
    M, L, t, M_acc, M_acc_max, M_acc_true = np.array([empty, empty, empty, empty, empty, empty])
    varis = np.full((n_t,xy[0],xy[1]), empty_vari)
    M[0] = (mc + me) * M_E
    L[0], varis[0] = L0, vari
    tn = 0
    t[0] = tn
    
    for n in range(n_t-1):
        
        Ln = L[n] * L_Sun
        #M[n+1] = M[n] + (Ln*R*dt)/(G*M[n])
        #M_acc[n] = (M[n+1]-M[n])/dt
        M_acc[n] = (Ln*R)/(G*M[n])
        
        if L09==True:
            Mdot_max, Sigma, R_sig = MARC.Max_Acc_Rate(Sigma_in,Ms,M[n],Td,Rp,r_in,r_out,alpha,num_sig,L09=True)
        else:
            Mdot_max, Sigma, R_sig = MARC.Max_Acc_Rate(Sigma_in,Ms,M[n],Td,Rp,r_in,r_out,alpha,num_sig)
        
        M_acc_max[n] = Mdot_max
        factor = 365*24*60*60/M_E
        print(M[n]/M_E,M_acc[n]*factor,M_acc_max[n]*factor)
        
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
        
        if M_acc[n] > M_acc_max[n]:
            print("Max Reached!")
            M_acc_true[n] = M_acc_max[n]
            M[n+1] = M[n] + (dt * M_acc_max[n])
        else:
            M_acc_true[n] = M_acc[n]
            M[n+1] = M[n] + (dt * M_acc[n])
        
        me += (M[n+1] - M[n]) / M_E
        
        m_grid *= M[n+1]/M[n]
        
        
        L[n+1], varis[n+1] = LFIO.L_Finder(mc,me,rp,pd,kd,Td,L[n],acc_l,n_i,m_grid, Plot=False)[:2]
        
        tn += dt
        t[n+1] = tn
    
    t /= (365*24*60*60)             # Convert time in s into years for plotting
    
    if M_acc[-1] > M_acc_max[-1]:
        print("Max Reached!")
        M_acc_true[-1] = M_acc_max[-1]
        M[-1] = M[-2] + (dt * M_acc_max[-1])
    else:
        M_acc_true[-1] = M_acc[-1]
        M[-1] = M[-2] + (dt * M_acc[-1])
    
    Mdots = [M_acc, M_acc_max, M_acc_true]
    
    return L, M, t, varis, Mdots
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 02:18:17 2023

@author: ollie
"""

import numpy as np
import scipy.constants as spc
import scipy.optimize as so
import matplotlib.pyplot as plt
import Data_Loader as DL





def fit_func(r,a,b):
    return a/np.sqrt(r) + b

def straight_func(r,b):
    return -1/2*r + b



def Sigma_Calc(Sigma_outer, q, H_r, nu, num_sig, Plot=False, Crida_Comparison=False):
    
    G = 1                           # Gravitational constant
    AU = spc.au * 1e2               # Astronomical unit [cm]
    
    Rp = 1
    Ms = 1
    Omega_p = np.sqrt(G*Ms/(Rp**3))
    R_H = Rp * (q/3)**(1/3)
    
    r = np.linspace(0.1, 3, int(num_sig))
    Sigma = np.zeros((int(num_sig)))
    
    Sigma[-1] = Sigma_outer
   
    for i in range(1, len(r)):
        delta_i = r[-i] - Rp
        delta_i_sign = delta_i/abs(delta_i)
        Omega_i = np.sqrt(G*Ms/(r[-i]**3))
        a_dd_i = ((1/8)*(abs(delta_i/R_H)**(-1.2))) + (200*(abs(delta_i/R_H)**(-10)))
        
        tg_i = 0.35*(q**2)*(Rp**5)*(Omega_p**2)*delta_i_sign*r[-i]*((1/delta_i)**4)
        t_f_nu_Omega_i = 3*nu*Omega_i/4
        H_r_term_i = (H_r**2)*R_H*Rp*(Omega_p**2)*r[-i]*a_dd_i
        t_f_nu_Omega_r_i = 3*nu*Omega_i*r[-i]/2
        
        Sigma_log = np.log(Sigma[-i]) + (((tg_i-t_f_nu_Omega_i)/(H_r_term_i+t_f_nu_Omega_r_i))*(r[-(i+1)]-r[-i]))
        Sigma[-(i+1)] = np.e ** Sigma_log
    
    
    if Plot == True:
        Figsize, lw, Ticksize, Labelsize, Legendfontsize = [7, 2, 20, 20, 15]
        
        plt.rc("figure", figsize=(Figsize,Figsize))
        plt.rc("xtick", labelsize=Ticksize)
        plt.rc("ytick", labelsize=Ticksize)
        plt.rc("axes", labelsize=Labelsize)
        plt.rc("legend", fontsize=Legendfontsize)
        
        plt.figure()
        plt.plot(r, Sigma, linewidth=lw, label="Analytical Model")
        
        if Crida_Comparison == True:
            Crida_Sigma_Data = DL.Data_Extractor_Sigma("Crida Data/Sigma (0.05, -5.5).csv")
            plt.plot(Crida_Sigma_Data[0], Crida_Sigma_Data[1], linewidth=lw, label="Crida Data")
        
        plt.xlabel("Radius / $R_p$")
        plt.ylabel("$\Sigma$")
        plt.legend()
        plt.show()
        plt.legend()
        ax = plt.gca()
        ax.set_box_aspect(1)
        plt.tight_layout()
        plt.show()
    
    r *= 5.2*AU
    Sigma *= 800
    
    return Sigma, r





def Max_Acc_Rate(Sigma_outer, Ms, Mp, Tp, Rp, r_in, r_out, alpha, num_sig, Const_Sigma=False, Plot=False, Crida_Comparison=False, L09=False):
    
    G = spc.G * 1e3                 # Gravitational constant
    kB = spc.k * 1e7                # Boltzmann constant
    mu = 2.3                        # Relative mass
    R = 8.31e7                      # Gas constant [erg/mol/K]
    R /= mu
    mp = spc.m_p * 1e3              # Proton mass [g]
    M_S = 1.989e33                  # Solar mass [g]
    M_E = 5.972e27                  # Earth mass [g]
    R_J = 6.9911e9
    AU = spc.au * 1e2               # Astronomical unit [cm]
    
    Ms *= M_S
    Rp *= AU
    
    if L09==True:
        cs2 = R * Tp           # Sound speed [cm/s]
        R_H = (5.2*AU) * (Mp/(3*M_S))**(1/3)
        r_in = 1e-10
        r_out = (G*Mp/((cs2/1.)+((G*Mp)/(0.25*R_H)))) / R_H
    
    
    if Crida_Comparison == True:
        H_r = 0.05
        nu = 3e14
        Mp = 1e-3 * Ms
        Sigma_outer = 1/np.sqrt(3)
        print("H/r =", H_r, "and \u03BD =", nu/1e14)
    
    else:
        cs_p = np.sqrt(R * Tp)
        H_p = cs_p*np.sqrt((Rp**3)/(G*Ms))
        H_r = H_p / Rp
        nu = alpha * cs_p * H_p
        #print("H/r =", H_r, "and \u03BD =", nu/1e14)
    
    Omega_p = np.sqrt(G*Ms/(Rp**3))
    R_H = Rp * (Mp/(3*Ms))**(1/3)
    
    
    if Mp == 0.:
        r_d = np.linspace(0.5*AU, 3*Rp, int(num_sig))
    
    else:
        r_d1 = np.linspace(0.5*AU, Rp-(r_out*R_H), int(num_sig/5)+1)[:-1]
        r_d2 = np.linspace(Rp-(r_out*R_H), Rp-(r_in*R_H), int(num_sig/5))
        r_d3 = np.linspace(Rp-(r_in*R_H), Rp+(r_in*R_H), int(num_sig/5)+2)[1:-1]
        r_d4 = np.linspace(Rp+(r_in*R_H), Rp+(r_out*R_H), int(num_sig/5))
        r_d5 = np.linspace(Rp+(r_out*R_H), 3*Rp, int(num_sig/5)+1)[1:]
        r_d = np.concatenate((r_d1,r_d2,r_d3,r_d4,r_d5))
    
    Sigma = np.zeros((len(r_d)))
    Mdot_maxs = []
    
    
    if Const_Sigma == True:
        Sigma[:] = Sigma_outer
        Mdot_max_ana_in = 2*Sigma_outer*Omega_p*(Rp**2)*(np.sqrt(1-(r_in*R_H/Rp))-np.sqrt(1-(r_out*R_H/Rp)))
        Mdot_max_ana_out = 2*Sigma_outer*Omega_p*(Rp**2)*(np.sqrt(1+(r_out*R_H/Rp))-np.sqrt(1+(r_in*R_H/Rp)))
        Mdot_max_ana = Mdot_max_ana_in + Mdot_max_ana_out
    
    else:
        Sigma[-1] = Sigma_outer
       
        for i in range(1, len(r_d)):
            delta_i = r_d[-i] - Rp
            delta_i_sign = delta_i/abs(delta_i)
            Omega_i = np.sqrt(G*Ms/(r_d[-i]**3))
            a_dd_i = ((1/8)*(abs(delta_i/R_H)**(-1.2))) + (200*(abs(delta_i/R_H)**(-10)))
            
            tg_i = 0.35*((Mp/Ms)**2)*(Rp**5)*(Omega_p**2)*delta_i_sign*r_d[-i]*((1/delta_i)**4)
            t_f_nu_Omega_i = 3*nu*Omega_i/4
            H_r_term_i = (H_r**2)*R_H*Rp*(Omega_p**2)*r_d[-i]*a_dd_i
            t_f_nu_Omega_r_i = 3*nu*Omega_i*r_d[-i]/2
            
            Sigma_log = np.log(Sigma[-i]) + (((tg_i-t_f_nu_Omega_i)/(H_r_term_i+t_f_nu_Omega_r_i))*(r_d[-(i+1)]-r_d[-i]))
            Sigma[-(i+1)] = np.e ** Sigma_log
    
    for i in range(int(num_sig/5)-1, int(2*num_sig/5)-1):
        Mdot_maxs.append(Sigma[i]*np.sqrt(G*Ms/r_d[i])*(r_d[i+1]-r_d[i]))
    for i in range(int(3*num_sig/5)-2, int(4*num_sig/5)-2):
        Mdot_maxs.append(Sigma[i]*np.sqrt(G*Ms/r_d[i])*(r_d[i+1]-r_d[i]))
    
    Mdot_max = np.sum(Mdot_maxs)
    
    
    if Plot == True:
        darkness = 0.5
        colour = "blue"
        Figsize, lw, Ticksize, Labelsize, Legendfontsize = [7, 2, 20, 20, 15]
        
        plt.rc("figure", figsize=(Figsize,Figsize))
        plt.rc("xtick", labelsize=Ticksize)
        plt.rc("ytick", labelsize=Ticksize)
        plt.rc("axes", labelsize=Labelsize)
        plt.rc("legend", fontsize=Legendfontsize)
        
        plt.figure()
        
        if Crida_Comparison == True:
            plt.plot(r_d/AU, Sigma, color="blue", linewidth=lw, label=f"0 Myr, ~{int(np.floor(Mp/M_E))} $M_{'E'}$")
            Crida_Sigma_Data = DL.Data_Extractor_Sigma("Crida Data/Sigma (0.05, -5.5).csv")
            plt.plot(Crida_Sigma_Data[0]*Rp/AU, Crida_Sigma_Data[1], color="orange", linewidth=lw, linestyle="--", label="Crida Data")
            
            #plt.axvspan(r_d[int(num_sig/5)]/AU, r_d[int(2*num_sig/5)-1]/AU, alpha=darkness, facecolor=colour)
            #plt.axvspan(r_d[int(3*num_sig/5)-1]/AU, r_d[int(4*num_sig/5)-1]/AU, alpha=darkness, facecolor=colour)
            #plt.axvline(r_d[int(num_sig/5)]/AU, color="blue", linewidth=lw, linestyle="--")
            #plt.axvline(r_d[int(2*num_sig/5)-1]/AU, color="blue", linewidth=lw, linestyle="--")
            #plt.axvline(r_d[int(3*num_sig/5)-1]/AU, color="blue", linewidth=lw, linestyle="--")
            #plt.axvline(r_d[int(4*num_sig/5)-1]/AU, color="blue", linewidth=lw, linestyle="--")
            plt.xlabel("Radius / AU")
            plt.ylabel("$\Sigma$ / $\Sigma_0$")
        
        else:
            if Mp == 0.:
                plt.loglog(r_d/AU, Sigma, color="blue", linewidth=lw, label=f"{Mp/M_E} $M_\oplus$")
                coeff = so.curve_fit(straight_func,np.log(r_d),np.log(Sigma))[0]
                Sigma_ana = straight_func(np.log(r_d),coeff[0])
                Sigma_ana = np.e**(Sigma_ana)
                plt.loglog(r_d/AU, Sigma_ana, color="orange", linewidth=lw, linestyle="--", label="$1/\sqrt{r}$ Fit")
            else:
                plt.plot(r_d/AU, Sigma, color="blue", linewidth=lw, label=f"0 Myr, {Mp/M_E} $M_\oplus$")
                plt.axvspan(r_d[int(num_sig/5)]/AU, r_d[int(2*num_sig/5)-1]/AU, alpha=darkness, facecolor=colour)
                plt.axvspan(r_d[int(3*num_sig/5)-1]/AU, r_d[int(4*num_sig/5)-1]/AU, alpha=darkness, facecolor=colour)
            plt.xlabel("Radius / AU")
            plt.ylabel("$\Sigma$ / g$\cdot$cm$^{-2}$")
        
        plt.legend()
        ax = plt.gca()
        ax.set_box_aspect(1)
        plt.tight_layout()
        plt.show()
    
    
    if Const_Sigma == True:
        return Mdot_max, Mdot_max_ana
    
    elif Crida_Comparison == True:
        return Sigma, r_d
    
    else:
        return Mdot_max, Sigma, r_d
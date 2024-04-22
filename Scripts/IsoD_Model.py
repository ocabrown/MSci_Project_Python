#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 19 16:52:33 2022

@author: ollie
"""

import numpy as np
import scipy.constants as spc
import scipy.optimize as so
import matplotlib.pyplot as plt





class Iso_Density_Giant_Planet:
    
    
    
    def __init__(self,mc,me,rp,pd,kd,Td,lp,n):
        
        # Constants:
        num_p = int(n)                  # Number of points along m
        self._n = num_p
        num_g = 1                       # Number of ghost cells on each side
        
        AU = spc.au * 1e2               # Astronomical unit [cm]
        G = spc.G * 1e3                 # Gravitational constant
        kB = spc.k * 1e7                # Boltzmann constant
        mu = 1.                        # Relative mass
        R = 8.31e7                      # Gas constant [erg/mol/K]
        R /= mu
        mp = spc.m_p * 1e3              # Proton mass [g]
        R_J = 6.9911e9                  # Jupiter radius [cm]
        L_Sun = 3.846e33                # Solar luminosity [erg/s]
        M_E = 5.972e27                  # Earth mass [g]
        M_S = 1.989e33
        
        cs2 = (7/5) * R * Td            # Sound speed [cm/s]
        
        self._lo = num_g
        self._hi = num_g + num_p - 1
        
        core_mass_max = mc * M_E        # Planet's core mass [M_E]
        env_mass_max = me * M_E         # Planet's envelope mass [M_E]
        Mp = core_mass_max + env_mass_max
        R_H = (5.2*AU) * (Mp/(3*M_S))**(1/3)
        rad_hi = (G*Mp/((cs2/1.)+((G*Mp)/(0.25*R_H))))
        #print(rad_hi/R_J)
        rad_hi = rp * R_J               # Planet radius [cm]
        den_hi = pd                     # Disc density [g/cm^3]
        opa_hi = kd                     # Disc opacity [cm^2/g]
        tem_hi = Td                     # Disc temperature [K]
        pre_hi = pd * R * Td            # Disc pressure [g/cm/s^2]
        lum = lp * L_Sun                # Luminosity [L_Sun]
        
        
        # Making empty arrays for other variables with boundary conditions
        empty = np.zeros((num_p + (2*num_g)))
        r, p, T, P, s = np.array([empty,empty,empty,empty,empty])
        r[self._hi] = rad_hi
        p[self._hi] = den_hi
        T[self._hi] = tem_hi
        P[self._hi] = pre_hi
        
        
        # Calculating envelope mass array so r has Uniform logarithmic spacing
        env_mass_min = 1e20         # [g]
        
        def m_r_func(r):
            integral = (4/3)*np.pi*(r**3)*den_hi
            C = env_mass_max - (4/3)*np.pi*(rad_hi**3)*den_hi
            return integral + C
        
        def m_min(r):
            return env_mass_min - m_r_func(r)
        
        
        rad_lo_guess = np.linspace(1e8, 1e11, 100000)
        m_mins = []
        
        for x in range(len(rad_lo_guess)):
            rad_lo_x = rad_lo_guess[x]
            m_min_x = m_min(rad_lo_x)
            m_mins.append(m_min_x)
        
        for x in range(len(m_mins)-1):
            if np.sign(m_mins[x]) != np.sign(m_mins[x+1]):
                rad_lo_lims = [rad_lo_guess[x], rad_lo_guess[x+1]]
                break
        
        rad_lo = so.brentq(m_min, rad_lo_lims[0], rad_lo_lims[1])
        
        r_lin_geom = np.geomspace(rad_lo, rad_hi, num=num_p)
        m = m_r_func(r_lin_geom)
        
        for i in range(num_g):
            m = np.insert(m, 0,      0.)
            m = np.insert(m, len(m), 0.)
        
        """
        
        # Uniform logarithmic spacing for me (envelope mass)
        env_mass_min = 1e20
        m1 = np.geomspace(env_mass_min, 0.9*env_mass_max, num=int(0.35*num_p))
        m2 = np.geomspace((0.9*env_mass_max), env_mass_max, num=int((0.65*num_p)+1))
        m2 = m2[1:]
        m = np.concatenate((m1,m2))
        m = np.geomspace(env_mass_min, env_mass_max, num=int(num_p))
        for i in range(num_g):
            m = np.insert(m, 0,      0.)
            m = np.insert(m, len(m), 0.)
        """
        
        self._const = [core_mass_max, rad_hi, opa_hi, tem_hi, lum, R_J, M_E]
        self._var_ini = [m, r, p, T, P, s]
    
    
    
    def Planet_Formation(self):
        
        # Constants:
        G = spc.G * 1e3                 # Gravitational constant
        a = 7.57e-15                    # erg/cm^3/K^4
        c = spc.c * 1e2                 # Light speed [cm/s]
        gamma = 7/5                     # Adiabatic constant
        M_core = self._const[0]
        r0 = self._const[1]
        k = self._const[2]
        T0 = self._const[3]
        l = self._const[4]
        
        var_new = self._var_ini.copy()
        m, r, p, T, P, s = var_new
        T_ana = T.copy()
        dP_dr_num = s.copy()
        dP_dr_ana = s.copy()
        
        
        for i in range(self._lo+1, self._hi+1):
            
            r[-(i+1)] = (r[-i]**3 + 3*(m[-(i+1)]-m[-i])/(4*np.pi*p[-i]))**(1/3)
            p[-(i+1)] = p[-i]
            T[-(i+1)] = (T[-i]**4 - 3*k*l*(m[-(i+1)]-m[-i])/(16*a*c*(np.pi**2)*r[-(i+1)]**4))**(1/4)
            P[-(i+1)] = P[-i] - G*(m[-(i+1)]+M_core)*(m[-(i+1)]-m[-i])/(4*np.pi*r[-(i+1)]**4)
            s[-(i+1)] = T[-(i+1)]/(P[-(i+1)]**((gamma-1)/gamma))
            
            T_ana[-(i+1)] = (T0**4 + ((3*k*l*p[-(i+1)]/(4*np.pi*a*c))*(1/r[-(i)] - 1/r0)))**(1/4)
            dP_dr_num[-(i+1)] = (P[-(i+1)]-P[-i])/(r[-(i+1)]-r[-i])
            dP_dr_ana[-(i+1)] = -G*(m[-(i+1)]+M_core)*p[-(i+1)]/r[-(i+1)]**2
        
        
        self._var = [m, r, p, T, P, s, T_ana, dP_dr_num, dP_dr_ana]
        
        return self._var
    
    
    
    def Plot(self, Variable, vs_m_or_r, Comparison_Data=None):
        
        R_J, M_E = self._const[5:]
        var_plot = self._var.copy()
        
        m, r, p, T, P, s, T_ana, dP_dr_num, dP_dr_ana = [var[self._lo:self._hi] for var in var_plot]
        
        # Numerical dP/dr straight line:
        min_pos = np.argmin(dP_dr_num)
        
        
        Figsize, lw, Ticksize, Labelsize, Legendfontsize = [7, 2, 20, 20, 12]
        
        plt.rc("figure", figsize=(Figsize,Figsize))
        plt.rc("xtick", labelsize=Ticksize)
        plt.rc("ytick", labelsize=Ticksize)
        plt.rc("axes", labelsize=Labelsize)
        plt.rc("legend", fontsize=Legendfontsize)
        
        plt.figure()
        
        # Plots against m
        if vs_m_or_r == "m":
            plt.xlabel("Mass / $M_{\oplus}$")
            
            if Variable == "r":
                plt.loglog(m/M_E, r/R_J, linewidth=lw)
                plt.ylabel("Radius / $R_{J}$")
                #plt.title("Radius vs. Mass")
            
            elif Variable == "p":
                plt.plot(m/M_E, p, linewidth=lw)
                plt.ylabel("Density / g$\cdot$cm$^{-3}$")
                #plt.title("Density vs. Mass")
            
            elif Variable == "T":
                plt.plot(m/M_E, T, linewidth=lw)
                plt.ylabel("Temperature / K")
                #plt.title("Temperature vs. Mass")
    
            elif Variable == "P":
                plt.plot(m/M_E, P, linewidth=lw)
                plt.ylabel("Pressure / dyne$\cdot$cm$^{-2}$")
                #plt.title("Pressure vs. Mass")
            
            elif Variable == "s":
                plt.plot(m/M_E, s, linewidth=lw)
                plt.ylabel(r"Measure of Entropy ($T/P^{(\gamma-1)/\gamma}$)")
                #plt.title("Measure of Entropy vs. Mass")
        
        # Plots against r
        elif vs_m_or_r == "r":
            plt.xlabel("Radius / $R_{J}$")
            plt.xticks(ticks=[])
            
            if Variable == "m":
                plt.loglog(r/R_J, m/M_E, linewidth=lw)
                plt.ylabel("Mass / $M_{\oplus}$")
                #plt.title("Mass vs. Radius")
            
            elif Variable == "p":
                plt.loglog(r/R_J, p, linewidth=lw)
                plt.ylabel("Density / g$\cdot$cm$^{-3}$")
                #plt.title("Density vs. Radius")
                plt.legend()
            
            elif Variable == "T":
                plt.plot(r/R_J, T, linewidth=lw)
                plt.xscale("log")
                plt.ylabel("Temperature / K")
                #plt.title("Temperature vs. Radius")
                plt.legend()
                
            elif Variable == "P":
                plt.plot(r/R_J, P, linewidth=lw)
                plt.xscale("log")
                plt.ylabel("Pressure / dyne$\cdot$cm$^{-2}$")
                #plt.title("Pressure vs. Radius")
                
            elif Variable == "s":
                plt.plot(r/R_J, s, linewidth=lw)
                plt.ylabel(r"Measure of Entropy ($T/P^{(\gamma-1)/\gamma}$)")
                #plt.title("Measure of Entropy vs. Radius")
            
            elif Variable == "T_comp":
                plt.plot(r/R_J, T, linewidth=lw, label="Numerical Solution")
                plt.plot(r/R_J, T_ana, "--", linewidth=lw, label="Analytical Solution")
                plt.ylabel("Temperature / K")
                plt.xscale("log")
                #plt.title("Temperature vs. Radius Comparison")
                plt.legend()
            
            elif Variable == "dP/dr":
                plt.plot(r/R_J, dP_dr_num/1e-9, linewidth=lw, label="Numerical Solution")
                plt.plot(r/R_J, dP_dr_ana/1e-9, "--", linewidth=lw, label="Analytical Solution")
                plt.ylabel("Pressure Gradient / ndyne$\cdot$cm$^{-3}$")
                plt.xscale("log")
                #plt.title("Pressure Gradient vs. Radius Comparison")
                plt.legend()
                
            elif Variable == "dP/dr_lower":
                plt.plot(r[:min_pos]/R_J, dP_dr_num[:min_pos]/1e-9, linewidth=lw, label="Numerical Solution")
                plt.plot(r[:min_pos]/R_J, dP_dr_ana[:min_pos]/1e-9, "--", linewidth=lw, label="Analytical Solution")
                plt.ylabel("Pressure Gradient / ndyne$\cdot$cm$^{-3}$")
                plt.xscale("log")
                #plt.title("Pressure Gradient vs. Radius - Low")
                plt.legend()
            
            elif Variable == "dP/dr_upper":
                plt.plot(r[min_pos:]/R_J, dP_dr_num[min_pos:]/1e-9, linewidth=lw, label="Numerical Solution")
                plt.plot(r[min_pos:]/R_J, dP_dr_ana[min_pos:]/1e-9, "--", linewidth=lw, label="Analytical Solution")
                plt.ylabel("Pressure Gradient / ndyne$\cdot$cm$^{-3}$")
                plt.xscale("log")
                #plt.title("Pressure Gradient vs. Radius - High")
                plt.legend()
            
            elif Variable == "sharedx":
                fig, ax1 = plt.subplots()
                color1 = 'red'
                color2 = 'blue'
                ax1.set_xlabel('Radius / $R_{J}$')
                ax1.set_ylabel("Temperaure / K", color=color1)
                l1 = ax1.plot(r/R_J, T, color=color1, linewidth=lw, label="Numerical Temperature Solution")
                ax1.plot(r/R_J, T_ana, ":", color=color2, linewidth=lw)
                ax1.tick_params(axis='y', labelcolor=color1)
                
                ax2 = ax1.twinx()
                
                color1 = 'orange'
                ax2.set_ylabel('Pressure Gradient / ndyne$\cdot$cm$^{-3}$', color=color1)
                l2 = ax2.plot(r/R_J, dP_dr_num/1e-9, color=color1, linewidth=lw, label="Numerical dP/dr Solution")
                l3 = ax2.plot(r/R_J, dP_dr_ana/1e-9, ":", color=color2, linewidth=lw, label="Analytical Solution")
                ax2.tick_params(axis='y', labelcolor=color1)
                fig.tight_layout()
                
                leg = l1 + l2 + l3
                labs = [l.get_label() for l in leg]
                plt.legend(leg, labs, loc=8)
        
            #plt.xticks([1.6e11, 1.9e11, 2.2e11])
        ax = plt.gca()
        ax.set_box_aspect(1)
        ax = plt.gca()
        ax.set_xticks([18,20,22,24,26,28,30])
        ax.set_xticklabels([18,20,22,24,26,28,30])
        plt.tight_layout()
        plt.show()
    
    
    
    def L2norm(self):
        
        var_conv = self._var.copy()
        m, r, p, T, P, s, T_ana, dP_dr_num, dP_dr_ana = [var[int(self._lo+self._n/10):int(self._hi-self._n/10)] for var in var_conv]
        
        min_pos = np.argmin(dP_dr_num)
        min_val = dP_dr_num[min_pos]
        
        L2norm_T = np.sqrt(sum(((T_ana - T)/T_ana)**2)/len(T))
        L2norms_T = ((T_ana - T)/T_ana)**2
        L2norm_dP_dr = np.sqrt(sum(((dP_dr_ana - dP_dr_num)/dP_dr_ana)**2)/len(dP_dr_num))
        #L2norm_dP_dr_lower = np.sqrt(sum(((dP_dr_ana[:min_pos] - dP_dr_num[:min_pos])/dP_dr_ana[:min_pos])**2)/len(dP_dr_num[:min_pos]))
        L2norm_dP_dr_lower = 1
        L2norm_dP_dr_upper = np.sqrt(sum(((dP_dr_ana[min_pos:] - dP_dr_num[min_pos:])/dP_dr_ana[min_pos:])**2)/len(dP_dr_num[min_pos:]))
        
        
        return L2norm_T, L2norm_dP_dr, L2norm_dP_dr_lower, L2norm_dP_dr_upper, min_val, L2norms_T
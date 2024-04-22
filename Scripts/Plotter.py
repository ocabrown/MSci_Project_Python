#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 20:30:24 2023

@author: ollie
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as spc
import Tau_Line_Finder as TLF



def Plot(Simulation_Data, y_variable, x_variable, Comparison_Data=None, Params=None, calc=False, Movshovitz=False, Evo=False, alpha_comp=False, acc_reg=False):
    
    R_J, M_E = 6.9911e9, 5.972e27       # Jupiter radius and Earth mass in cgs
    L_J, L_S = 5e24, 3.846e33     # Jupiter and Solar luminosity [erg/s]
    
    
    plt.figure()
    
    if x_variable == "calc":
        
        L, M, t, varis, Mdots = Simulation_Data[0]
        L_SB = []
        r_tau = []
        T_tau = []
        
        for n in range(len(varis)):
            lo = 1
            hi = -2
            
            r_plot = varis[n][1][lo:hi]
            p_plot = varis[n][2][lo:hi]
            k_plot = varis[n][3][lo:hi]
            T_plot = varis[n][4][lo:hi]
            
            tau = 1
            tau_i = TLF.Tau_Line_Finder(r_plot, p_plot, k_plot, tau)
            
            r_tau_i = r_plot[-tau_i]
            r_tau.append(r_tau_i)
            T_tau_i = T_plot[-tau_i]
            T_tau.append(T_tau_i)
            
            import scipy.constants as spc
            sigma_SB = spc.Stefan_Boltzmann * 1e3
            
            L_SB_n = 4 * np.pi * (r_tau_i**2) * sigma_SB * (T_tau_i**4) / L_S
            L_SB.append(L_SB_n)
            
        L_SB = np.array(L_SB)
        r_tau = np.array(r_tau)
        T_tau = np.array(T_tau)
        
        x = [r_tau,T_tau,L_SB]
        return x
    
    if x_variable == "t":
        
        Figsize, lw, Ticksize, Labelsize, Legendfontsize, Titlesize = Params
        
        plt.rc("figure", figsize=(Figsize,Figsize))
        plt.xlabel("Radius / $R_{J}$")
        plt.rc("xtick", labelsize=Ticksize)
        plt.rc("ytick", labelsize=Ticksize)
        plt.rc("axes", labelsize=Labelsize)
        plt.rc("legend", fontsize=Legendfontsize)
        
        L, M, t, varis, Mdots = Simulation_Data[0]
        factor = 365*24*60*60/M_E
        Mdots = np.array(Mdots)
        M_acc, M_acc_max, M_acc_true = Mdots*factor
        
        if calc != False:
            r_tau, T_tau, L_SB = calc[0], calc[1], calc[2]
            import scipy.constants as spc
            sigma_SB = spc.Stefan_Boltzmann * 1e3
            L_SB1 = 4 * np.pi * (r_tau**2) * sigma_SB * (T_tau**4) / L_S
            
        else:
            L_SB = []
            r_tau = []
            T_tau = []
            
            for n in range(len(varis)):
                lo = 1
                hi = -2
                
                r_plot = varis[n][1][lo:hi]
                p_plot = varis[n][2][lo:hi]
                k_plot = varis[n][3][lo:hi]
                T_plot = varis[n][4][lo:hi]
                
                tau = 1
                tau_i = TLF.Tau_Line_Finder(r_plot, p_plot, k_plot, tau)
                
                r_tau_i = r_plot[-tau_i]
                r_tau.append(r_tau_i)
                T_tau_i = T_plot[-tau_i]
                T_tau.append(T_tau_i)
                
                import scipy.constants as spc
                sigma_SB = spc.Stefan_Boltzmann * 1e3
                
                L_SB_n = 4 * np.pi * (r_tau_i**2) * sigma_SB * (T_tau_i**4) / L_S
                L_SB.append(L_SB_n)
                
            L_SB = np.array(L_SB)
            r_tau = np.array(r_tau)
            T_tau = np.array(T_tau)
        
        
        if y_variable == "L":
            plt.plot(t/1e6, L/1e-5, linewidth=lw)
            plt.ylabel("Luminosity / 10$^{-5}$ $L_{\odot}$")
            #plt.ylabel("Luminosity / $L_{\odot}$")
            #plt.title("Luminosity vs. Time")
        
        if y_variable == "L_t":
            plt.plot(t/1e6, L_SB/(1e-5), linewidth=lw, color="orange", linestyle="-", label="Stefan-Boltzmann Law")
            #plt.plot(t/1e6, L/1e-5*0.55, linewidth=lw, color="red", linestyle="--", label="0.55x Numerical")
            plt.plot(t/1e6, L/1e-5, linewidth=lw, color="red", linestyle="-", label="Numerical")
            plt.plot(t/1e6, (L_SB - L)/1e-5, linewidth=lw, color="blue", linestyle="-")
            plt.ylabel("Luminosity / 10$^{-5}$ $L_{\odot}$")
            plt.legend()
            #plt.title("Stefan-Boltzmann Luminosity vs. Radius")
        
        if y_variable == "L_comp":
            plt.plot(t/1e6, (L_SB - L)/1e-5, linewidth=lw, color="blue", linestyle="-")
            plt.ylabel("Luminosity / 10$^{-5}$ $L_{\odot}$")
            plt.legend()
        
        if y_variable == "T_t,r_t":
            Figsize, lw, Ticksize, Labelsize, Legendfontsize, Titlesize = Params
            
            plt.rc("figure", figsize=(Figsize,Figsize))
            plt.rc("xtick", labelsize=Ticksize)
            plt.rc("ytick", labelsize=Ticksize)
            plt.rc("axes", labelsize=Labelsize)
            plt.rc("legend", fontsize=Legendfontsize)
            
            fig, ax1 = plt.subplots()
            color1 = 'red'
            color2 = 'orange'
            color3 = 'blue'
            
            ax1.set_xlabel('Time / Myrs')
            ax1.set_ylabel(r'$r_{\tau=1}$ / $R_J$', color=color1)
            ax1.plot(t/1e6, r_tau / R_J, linewidth=lw, color=color1)
            ax1.tick_params(axis='y')
            
            ax2 = ax1.twinx()
            ax2.set_ylabel(r'$T_{\tau=1}$ / K', color=color2)
            ax2.plot(t/1e6, T_tau, color=color2, linewidth=lw)
            ax2.tick_params(axis='y')
            fig.tight_layout()
            
            ax3 = ax1.twinx()
            ax3.set_ylabel("Luminosity / 10$^{-5}$ $L_{\odot}$", color=color3)
            ax3.spines.right.set_position(("axes", 1.3))
            ax3.plot(t/1e6, L/1e-5, color=color3, linewidth=lw, linestyle="--")
            ax3.tick_params(axis='y')
            fig.tight_layout()
        
        if y_variable == "T_t4,r_t2":
            Figsize, lw, Ticksize, Labelsize, Legendfontsize, Titlesize = Params
            
            plt.rc("figure", figsize=(Figsize,Figsize))
            plt.rc("xtick", labelsize=Ticksize)
            plt.rc("ytick", labelsize=Ticksize)
            plt.rc("axes", labelsize=Labelsize)
            plt.rc("legend", fontsize=Legendfontsize)
            
            fig, ax1 = plt.subplots()
            color1 = 'red'
            color2 = 'orange'
            color3 = 'blue'
            
            ax1.set_xlabel('Time / Myrs')
            ax1.set_ylabel(r'$r_{\tau=1}^2$ / $R_J^2$', color=color1)
            ax1.plot(t/1e6, (r_tau/R_J)**2, linewidth=lw, color=color1)
            ax1.tick_params(axis='y')
            
            ax2 = ax1.twinx()
            ax2.set_ylabel(r'$T_{\tau=1}^4$ / $10^{10}$ K$^4$', color=color2)
            ax2.plot(t/1e6, T_tau**4/1e10, color=color2, linewidth=lw)
            ax2.tick_params(axis='y')
            fig.tight_layout()
            
            ax3 = ax1.twinx()
            ax3.set_ylabel(r'$r_{\tau=1}^2 \cdot T_{\tau=1}^4$ / $10^{32}$ cm$^2 \cdot$K$^4$', color=color3)
            ax3.spines.right.set_position(("axes", 1.2))
            ax3.plot(t/1e6, ((r_tau**2) * (T_tau**4))/1e32, color=color3, linewidth=lw, linestyle="--")
            ax3.tick_params(axis='y')
            fig.tight_layout()
        
        
        if y_variable == "M":
            plt.plot(t/1e6, M/M_E, linewidth=lw)
            plt.ylabel("Total Mass / $M_{\oplus}$")
            #plt.title("Total Mass vs. Time")
        
        if y_variable == "Mdot":
            plt.plot(t[:-1]/1e6, M_acc_true[:-1], linestyle="-", linewidth=lw, color="blue", label="$\dot{M}$")
            plt.plot(t[:-1]/1e6, M_acc[:-1], linestyle="--", linewidth=lw, color="orange", label="$\dot{M}_{basic}$")
            plt.plot(t[:-1]/1e6, M_acc_max[:-1], linestyle="--", linewidth=lw, color="red", label="$\dot{M}_{lim}$")
            
            plt.yscale("log")
            plt.ylabel("Accretion Rate / $M_{\oplus}\cdot$yr$^{-1}$")
            plt.legend()
            #plt.title("Accretion Rate vs. Time")
        
        if y_variable == "Mdot+M":
            Figsize, lw, Ticksize, Labelsize, Legendfontsize, Titlesize = Params
            
            plt.rc("figure", figsize=(Figsize,Figsize))
            plt.rc("xtick", labelsize=Ticksize)
            plt.rc("ytick", labelsize=Ticksize)
            plt.rc("axes", labelsize=Labelsize)
            plt.rc("legend", fontsize=Legendfontsize)
            
            fig, ax1 = plt.subplots()
            color1 = 'red'
            color2 = 'orange'
            
            ax1.set_xlabel('Time / Myrs')
            ax1.set_ylabel('Accretion Rate / $M_{\oplus}\cdot$yr$^{-1}$', color=color1)
            l1 = ax1.plot(t[:-1]/1e6, M_acc[:-1], linestyle=":", linewidth=lw, color=color1, label="$\dot{M}_{basic}$")
            l2 = ax1.plot(t[:-1]/1e6, M_acc_max[:-1], linestyle="--", linewidth=lw, color=color1, label="$\dot{M}_{lim}$")
            l3 = ax1.plot(t[:-1]/1e6, M_acc_true[:-1], linestyle="-", linewidth=lw, color=color1, label="$\dot{M}$")
            plt.yscale("log")
            ax1.tick_params(axis='y')
            
            ax2 = ax1.twinx()
            ax2.set_ylabel('Total Mass / $M_{\oplus}$', color=color2)
            ax2.plot(t[:]/1e6, M[:]/M_E, color=color2, linewidth=lw)
            ax2.tick_params(axis='y')
            fig.tight_layout()
            
            leg = l1 + l2 + l3
            labs = [l.get_label() for l in leg]
            plt.legend(leg, labs, loc=7)
        
        if y_variable == "L+M":
            Figsize, lw, Ticksize, Labelsize, Legendfontsize, Titlesize = Params
            
            plt.rc("figure", figsize=(Figsize,Figsize))
            plt.rc("xtick", labelsize=Ticksize)
            plt.rc("ytick", labelsize=Ticksize)
            plt.rc("axes", labelsize=Labelsize)
            plt.rc("legend", fontsize=Legendfontsize)
            
            fig, ax1 = plt.subplots()
            color1 = 'red'
            color2 = 'orange'
            
            ax1.set_xlabel('Time / Myrs')
            ax1.set_ylabel('Luminosity / 10$^{-5}$ $L_{\odot}$', color=color1)
            #ax1.set_ylabel('Luminosity / $L_{\odot}$', color=color1)
            ax1.plot(t[:]/1e6, L[:]/1e-5, color=color1, linewidth=lw)
            ax1.tick_params(axis='y')
            
            ax2 = ax1.twinx()
            ax2.set_ylabel('Total Mass / $M_{\oplus}$', color=color2)
            ax2.plot(t[:]/1e6, M[:]/M_E, color=color2, linewidth=lw)
            ax2.tick_params(axis='y')
            fig.tight_layout()
        
        plt.xlabel("Time / Myrs")
        #plt.xscale("log")
        ax = plt.gca()
        ax.set_box_aspect(1)
        ax = plt.gca()
        #ax.set_xticks([0.01,0.1,1,2,3])
        #ax.set_xticklabels([0.01,0.1,1,2,3])
        plt.tight_layout()
        fig.tight_layout()
        plt.show()
    
    
    
    elif len(Simulation_Data) == 4 and x_variable != "t":
        var_plot0 = Simulation_Data[0]._var.copy()
        var_plot1 = Simulation_Data[1]._var.copy()
        var_plot2 = Simulation_Data[2]._var.copy()
        var_plot3 = Simulation_Data[3]._var.copy()
        m0, r0, p0, k0, T0, P0, s0 = [var[Simulation_Data[0]._lo:Simulation_Data[0]._hi] for var in var_plot0]
        m1, r1, p1, k1, T1, P1, s1 = [var[Simulation_Data[1]._lo:Simulation_Data[1]._hi] for var in var_plot1]
        m2, r2, p2, k2, T2, P2, s2 = [var[Simulation_Data[2]._lo:Simulation_Data[2]._hi] for var in var_plot2]
        m3, r3, p3, k3, T3, P3, s3 = [var[Simulation_Data[3]._lo:Simulation_Data[3]._hi] for var in var_plot3]
        
        Figsize, lw, Ticksize, Labelsize, Legendfontsize, Titlesize = Params
        
        
        plt.rc("figure", figsize=(Figsize,Figsize))
        plt.xlabel("Radius / $R_{J}$")
        plt.rc("xtick", labelsize=Ticksize)
        plt.rc("ytick", labelsize=Ticksize)
        plt.rc("axes", labelsize=Labelsize)
        plt.rc("legend", fontsize=Legendfontsize)
        
        
        if y_variable != "p" and y_variable != "k" and y_variable != "T":
            print("Haven't written code to plot", y_variable, "for the report!")
        
        
        elif y_variable == "p":
            
            for x in range(len(Simulation_Data[0]._boundary_i)-1):
                high = -Simulation_Data[0]._boundary_i[x+1][0]
                low = -Simulation_Data[0]._boundary_i[x][0]
                if Simulation_Data[0]._boundary_i[x][1] == 1:
                    if x == 0 or x == 1:
                        plt.loglog(r0[high:low]/R_J, p0[high:low], linewidth=lw[0], linestyle="-", color="royalblue", label="Constant Opacity = 10$^{-4}$ cm$^{2}$g$^{-1}$")
                    else:
                        plt.loglog(r0[high:low]/R_J, p0[high:low], linewidth=lw[0], linestyle="-", color="royalblue")
                if Simulation_Data[0]._boundary_i[x][1] == 0:
                    plt.loglog(r0[high:low]/R_J, p0[high:low], linewidth=lw[1], linestyle="-", color="royalblue")
            
            for x in range(len(Simulation_Data[1]._boundary_i)-1):
                high = -Simulation_Data[1]._boundary_i[x+1][0]
                low = -Simulation_Data[1]._boundary_i[x][0]
                if Simulation_Data[1]._boundary_i[x][1] == 1:
                    if x == 0 or x == 1:
                        plt.loglog(r1[high:low]/R_J, p1[high:low], linewidth=lw[0], linestyle="-", color="forestgreen", label="Movshovitz Opacity")
                    else:
                        plt.loglog(r1[high:low]/R_J, p1[high:low], linewidth=lw[0], linestyle="-", color="forestgreen")
                if Simulation_Data[1]._boundary_i[x][1] == 0:
                    plt.loglog(r1[high:low]/R_J, p1[high:low], linewidth=lw[1], linestyle="-", color="forestgreen")
            
            for x in range(len(Simulation_Data[2]._boundary_i)-1):
                high = -Simulation_Data[2]._boundary_i[x+1][0]
                low = -Simulation_Data[2]._boundary_i[x][0]
                if Simulation_Data[2]._boundary_i[x][1] == 1:
                    if x == 0 or x == 1:
                        plt.loglog(r2[high:low]/R_J, p2[high:low], linewidth=lw[0], linestyle="--", color="firebrick", label="Freedman Opacity (L = 10L$_{0}$)")
                    else:
                        plt.loglog(r2[high:low]/R_J, p2[high:low], linewidth=lw[0], linestyle="--", color="firebrick")
                if Simulation_Data[2]._boundary_i[x][1] == 0:
                    plt.loglog(r2[high:low]/R_J, p2[high:low], linewidth=lw[1], linestyle="--", color="firebrick")
            
            for x in range(len(Simulation_Data[3]._boundary_i)-1):
                high = -Simulation_Data[3]._boundary_i[x+1][0]
                low = -Simulation_Data[3]._boundary_i[x][0]
                if Simulation_Data[3]._boundary_i[x][1] == 1:
                    if x == 0 or x == 1:
                        plt.loglog(r3[high:low]/R_J, p3[high:low], linewidth=lw[0], linestyle="-", color="firebrick", label="Freedman Opacity (L = L$_{0}$)")
                    else:
                        plt.loglog(r3[high:low]/R_J, p3[high:low], linewidth=lw[0], linestyle="-", color="firebrick")
                if Simulation_Data[3]._boundary_i[x][1] == 0:
                    plt.loglog(r3[high:low]/R_J, p3[high:low], linewidth=lw[1], linestyle="-", color="firebrick")
            
            #for x in range(len(Simulation_Data[3]._boundary_i)-1):
            #    high = -Simulation_Data[3]._boundary_i[x+1][0]
            #    low = -Simulation_Data[3]._boundary_i[x][0]
            #    if Simulation_Data[3]._boundary_i[x][1] == 1:
            #        if x == 0 or x == 1:
            #            plt.loglog(r3[high:low]/R_J, p3[high:low], linewidth=lw[0], linestyle="-", color="darkblue", label="Constant Opacity = 10$^{-?}$ cm$^{2}$g$^{-1}$")
            #        else:
            #            plt.loglog(r3[high:low]/R_J, p3[high:low], linewidth=lw[0], linestyle="-", color="darkblue")
            #    if Simulation_Data[3]._boundary_i[x][1] == 0:
            #        plt.loglog(r3[high:low]/R_J, p3[high:low], linewidth=lw[1], linestyle="-", color="darkblue")
            
            plt.ylabel("Density / g$\cdot$cm$^{-3}$")
            #plt.title("Density vs. Radius")
        
        
        elif y_variable == "k":
            
            for x in range(len(Simulation_Data[0]._boundary_i)-1):
                high = -Simulation_Data[0]._boundary_i[x+1][0]
                low = -Simulation_Data[0]._boundary_i[x][0]
                if Simulation_Data[0]._boundary_i[x][1] == 1:
                    if x == 0 or x == 1:
                        plt.loglog(r0[high:low]/R_J, k0[high:low], linewidth=lw[0], linestyle="-", color="royalblue", label="Constant Opacity = 10$^{-4}$ cm$^{2}$g$^{-1}$")
                    else:
                        plt.loglog(r0[high:low]/R_J, k0[high:low], linewidth=lw[0], linestyle="-", color="royalblue")
                if Simulation_Data[0]._boundary_i[x][1] == 0:
                    plt.loglog(r0[high:low]/R_J, k0[high:low], linewidth=lw[1], linestyle="-", color="royalblue")
            
            for x in range(len(Simulation_Data[1]._boundary_i)-1):
                high = -Simulation_Data[1]._boundary_i[x+1][0]
                low = -Simulation_Data[1]._boundary_i[x][0]
                if Simulation_Data[1]._boundary_i[x][1] == 1:
                    if x == 0 or x == 1:
                        plt.loglog(r1[high:low]/R_J, k1[high:low], linewidth=lw[0], linestyle="-", color="forestgreen", label="Movshovitz Opacity")
                    else:
                        plt.loglog(r1[high:low]/R_J, k1[high:low], linewidth=lw[0], linestyle="-", color="forestgreen")
                if Simulation_Data[1]._boundary_i[x][1] == 0:
                    plt.loglog(r1[high:low]/R_J, k1[high:low], linewidth=lw[1], linestyle="-", color="forestgreen")
            
            for x in range(len(Simulation_Data[2]._boundary_i)-1):
                high = -Simulation_Data[2]._boundary_i[x+1][0]
                low = -Simulation_Data[2]._boundary_i[x][0]
                if Simulation_Data[2]._boundary_i[x][1] == 1:
                    if x == 0 or x == 1:
                        plt.loglog(r2[high:low]/R_J, k2[high:low], linewidth=lw[0], linestyle="--", color="firebrick", label="Freedman Opacity (L = 10L$_{0}$)")
                    else:
                        plt.loglog(r2[high:low]/R_J, k2[high:low], linewidth=lw[0], linestyle="--", color="firebrick")
                if Simulation_Data[2]._boundary_i[x][1] == 0:
                    plt.loglog(r2[high:low]/R_J, k2[high:low], linewidth=lw[1], linestyle="--", color="firebrick")
            
            for x in range(len(Simulation_Data[3]._boundary_i)-1):
                high = -Simulation_Data[3]._boundary_i[x+1][0]
                low = -Simulation_Data[3]._boundary_i[x][0]
                if Simulation_Data[3]._boundary_i[x][1] == 1:
                    if x == 0 or x == 1:
                        plt.loglog(r3[high:low]/R_J, k3[high:low], linewidth=lw[0], linestyle="-", color="firebrick", label="Freedman Opacity (L = L$_{0}$)")
                    else:
                        plt.loglog(r3[high:low]/R_J, k3[high:low], linewidth=lw[0], linestyle="-", color="firebrick")
                if Simulation_Data[3]._boundary_i[x][1] == 0:
                    plt.loglog(r3[high:low]/R_J, k3[high:low], linewidth=lw[1], linestyle="-", color="firebrick")
            
            #for x in range(len(Simulation_Data[3]._boundary_i)-1):
            #    high = -Simulation_Data[3]._boundary_i[x+1][0]
            #    low = -Simulation_Data[3]._boundary_i[x][0]
            #    if Simulation_Data[3]._boundary_i[x][1] == 1:
            #        if x == 0 or x == 1:
            #            plt.loglog(r3[high:low]/R_J, k3[high:low], linewidth=lw[0], linestyle="-", color="darkblue", label="Constant Opacity = 10$^{-?}$ cm$^{2}$g$^{-1}$")
            #        else:
            #            plt.loglog(r3[high:low]/R_J, k3[high:low], linewidth=lw[0], linestyle="-", color="darkblue")
            #    if Simulation_Data[3]._boundary_i[x][1] == 0:
            #        plt.loglog(r3[high:low]/R_J, k3[high:low], linewidth=lw[1], linestyle="-", color="darkblue")
            
            plt.ylabel("Opacity / cm$^{2}$$\cdot$g$^{-1}$")
            #plt.title("Opacity vs. Radius")
        
        
        elif y_variable == "T":
            
            for x in range(len(Simulation_Data[0]._boundary_i)-1):
                high = -Simulation_Data[0]._boundary_i[x+1][0]
                low = -Simulation_Data[0]._boundary_i[x][0]
                if Simulation_Data[0]._boundary_i[x][1] == 1:
                    if x == 0 or x == 1:
                        plt.loglog(r0[high:low]/R_J, T0[high:low], linewidth=lw[0], linestyle="-", color="royalblue", label="Constant Opacity = 10$^{-4}$ cm$^{2}$g$^{-1}$")
                    else:
                        plt.loglog(r0[high:low]/R_J, T0[high:low], linewidth=lw[0], linestyle="-", color="royalblue")
                if Simulation_Data[0]._boundary_i[x][1] == 0:
                    plt.loglog(r0[high:low]/R_J, T0[high:low], linewidth=lw[1], linestyle="-", color="royalblue")
            
            for x in range(len(Simulation_Data[1]._boundary_i)-1):
                high = -Simulation_Data[1]._boundary_i[x+1][0]
                low = -Simulation_Data[1]._boundary_i[x][0]
                if Simulation_Data[1]._boundary_i[x][1] == 1:
                    if x == 0 or x == 1:
                        plt.loglog(r1[high:low]/R_J, T1[high:low], linewidth=lw[0], linestyle="-", color="forestgreen", label="Movshovitz Opacity")
                    else:
                        plt.loglog(r1[high:low]/R_J, T1[high:low], linewidth=lw[0], linestyle="-", color="forestgreen")
                if Simulation_Data[1]._boundary_i[x][1] == 0:
                    plt.loglog(r1[high:low]/R_J, T1[high:low], linewidth=lw[1], linestyle="-", color="forestgreen")
            
            for x in range(len(Simulation_Data[2]._boundary_i)-1):
                high = -Simulation_Data[2]._boundary_i[x+1][0]
                low = -Simulation_Data[2]._boundary_i[x][0]
                if Simulation_Data[2]._boundary_i[x][1] == 1:
                    if x == 0 or x == 1:
                        plt.loglog(r2[high:low]/R_J, T2[high:low], linewidth=lw[0], linestyle="--", color="firebrick", label="Freedman Opacity (L = 10L$_{0}$)")
                    else:
                        plt.loglog(r2[high:low]/R_J, T2[high:low], linewidth=lw[0], linestyle="--", color="firebrick")
                if Simulation_Data[2]._boundary_i[x][1] == 0:
                    plt.loglog(r2[high:low]/R_J, T2[high:low], linewidth=lw[1], linestyle="--", color="firebrick")
            
            for x in range(len(Simulation_Data[3]._boundary_i)-1):
                high = -Simulation_Data[3]._boundary_i[x+1][0]
                low = -Simulation_Data[3]._boundary_i[x][0]
                if Simulation_Data[3]._boundary_i[x][1] == 1:
                    if x == 0 or x == 1:
                        plt.loglog(r3[high:low]/R_J, T3[high:low], linewidth=lw[0], linestyle="-", color="firebrick", label="Freedman Opacity (L = L$_{0}$)")
                    else:
                        plt.loglog(r3[high:low]/R_J, T3[high:low], linewidth=lw[0], linestyle="-", color="firebrick")
                if Simulation_Data[3]._boundary_i[x][1] == 0:
                    plt.loglog(r3[high:low]/R_J, T3[high:low], linewidth=lw[1], linestyle="-", color="firebrick")
            
            #for x in range(len(Simulation_Data[3]._boundary_i)-1):
            #    high = -Simulation_Data[3]._boundary_i[x+1][0]
            #    low = -Simulation_Data[3]._boundary_i[x][0]
            #    if Simulation_Data[3]._boundary_i[x][1] == 1:
            #        if x == 0 or x == 1:
            #            plt.loglog(r3[high:low]/R_J, T3[high:low], linewidth=lw[0], linestyle="-", color="darkblue", label="Constant Opacity = 10$^{-?}$ cm$^{2}$g$^{-1}$")
            #        else:
            #            plt.loglog(r3[high:low]/R_J, T3[high:low], linewidth=lw[0], linestyle="-", color="darkblue")
            #    if Simulation_Data[3]._boundary_i[x][1] == 0:
            #        plt.loglog(r3[high:low]/R_J, T3[high:low], linewidth=lw[1], linestyle="-", color="darkblue")
            
            plt.ylabel("Temperature / K")
            #plt.xscale("log")
            #plt.title("Temperature vs. Radius")
        
        plt.legend(loc="upper right")
        ax = plt.gca()
        ax.set_box_aspect(1)
        plt.tight_layout()
        plt.show()
    
    
    
    elif Evo == True:
        if len(Simulation_Data) == 3:
            
            Figsize, lw, Ticksize, Labelsize, Legendfontsize, Titlesize = Params
            
            plt.rc("figure", figsize=(Figsize,Figsize))
            plt.xlabel("Time / Myr")
            plt.rc("xtick", labelsize=Ticksize)
            plt.rc("ytick", labelsize=Ticksize)
            plt.rc("axes", labelsize=Labelsize)
            plt.rc("legend", fontsize=Legendfontsize)
            plt.rc("axes", titlesize=Titlesize)
            
            
            L0, M0, t0, varis0, Mdots0 = Simulation_Data[0]
            L1, M1, t1, varis1, Mdots1 = Simulation_Data[1]
            L2, M2, t2, varis2, Mdots2 = Simulation_Data[2]
            
            factor = 365*24*60*60/M_E
            Mdots0 = np.array(Mdots0)
            Mdots1 = np.array(Mdots1)
            Mdots2 = np.array(Mdots2)
            
            M_acc0, M_acc_max0, M_acc_true0 = Mdots0*factor
            M_acc1, M_acc_max1, M_acc_true1 = Mdots1*factor
            M_acc2, M_acc_max2, M_acc_true2 = Mdots2*factor
            
            #i_cut = int(np.where(t0==2.49e6)[0])
            #print(i_cut)
            #print(t0[i_cut])
            #klog0 = np.floor(np.log10(varis0[0][3][1])*10)/10
            #klog1 = np.floor(np.log10(varis1[0][3][1])*10)/10
            #klog2 = np.floor(np.log10(varis2[0][3][1])*10)/10
            
            #467
            if y_variable == "L":
                plt.plot(t0[:]/1e6, L0[:]/1e-5, linewidth=lw, color="blue", label="$\kappa$ = $\kappa_0$")
                plt.plot(t1[:]/1e6, L1[:]/1e-5, linewidth=lw, color="red", label=r"$\kappa$ = $\kappa_0$/6")
                plt.plot(t2[:]/1e6, L2[:]/1e-5, linewidth=lw, color="orange", label=r"$\kappa$ = $\kappa_0$/12")
                #plt.plot(t2[:467]/1e6, L2[:467]/1e-6, linewidth=lw, color="firebrick", label=f"Constant Opacity = $10^{{{klog2}}}/10^{{{klog1}}}\kappa_0$")
                plt.yscale("log")
                plt.ylabel("Luminosity / 10$^{-5}$ $L_{\odot}$")
                #plt.title("Luminosity vs. Time")
                ax = plt.gca()
                ax.set_yticks([1,2,3,4,5,6])
                ax.set_yticklabels([1,2,3,4,5,6])
            
            if y_variable == "M":
                plt.plot(t0[:]/1e6, M0[:]/M_E, linewidth=lw, color="blue", label="$\kappa$ = $\kappa_0$")
                plt.plot(t1[:]/1e6, M1[:]/M_E, linewidth=lw, color="red", label=r"$\kappa$ = $\kappa_0$/6")
                plt.plot(t2[:]/1e6, M2[:]/M_E, linewidth=lw, color="orange", label=r"$\kappa$ = $\kappa_0$/12")
                #plt.yscale("log")
                plt.ylabel("Total Mass / $M_{\oplus}$")
                #plt.title("Total Mass vs. Time")
            
            if y_variable == "Mdot":
                if acc_reg == True:
                    plt.plot(t0[:-1]/1e6, M_acc0[:-1], linestyle="-", linewidth=lw, color="grey", label=r"$\dot{M}_{basic}$")
                    plt.plot(t1[:-1]/1e6, M_acc_max1[:-1], linestyle="--", linewidth=lw, color="blue", label=r"$\dot{M}_{lim}$, $\alpha$ = 10$^{-2}$, T12")
                    plt.plot(t0[:-1]/1e6, M_acc_max0[:-1], linestyle="--", linewidth=lw, color="red", label=r"$\dot{M}_{lim}$, $\alpha$ = 10$^{-5}$, T12")
                    plt.plot(t2[:-1]/1e6, M_acc_max2[:-1], linestyle="--", linewidth=lw, color="orange", label=r"$\dot{M}_{lim}$, $\alpha$ = 10$^{-5}$, L09")
                
                elif alpha_comp == True:
                    #plt.plot(t0[:-1]/1e6, M_acc_true0[:-1], linestyle="-", linewidth=lw, color="blue", label=r"$\dot{M}$, $\alpha$ = 10$^{-5}$")
                    plt.plot(t0[:-1]/1e6, M_acc0[:-1], linestyle="-", linewidth=lw, color="blue", label=r"$\dot{M}_{basic}$, $\alpha$ = 10$^{-5}$")
                    plt.plot(t0[:-1]/1e6, M_acc_max0[:-1], linestyle="--", linewidth=lw, color="blue", label=r"$\dot{M}_{lim}$, $\alpha$ = 10$^{-5}$")
                    #plt.plot(t1[:-1]/1e6, M_acc_true1[:-1], linestyle="-", linewidth=lw, color="red", label=r"$\dot{M}$, $\alpha$ = 10$^{-4}$")
                    plt.plot(t1[:-1]/1e6, M_acc1[:-1], linestyle="-", linewidth=lw, color="red", labelr="$\dot{M}_{basic}$, $\alpha$ = 10$^{-4}$")
                    plt.plot(t1[:-1]/1e6, M_acc_max1[:-1], linestyle="--", linewidth=lw, color="red", label=r"$\dot{M}_{lim}$, $\alpha$ = 10$^{-4}$")
                    #plt.plot(t2[:-1]/1e6, M_acc_true2[:-1], linestyle="-", linewidth=lw, color="orange", label=r"$\dot{M}$, $\alpha$ = 10$^{-3}$")
                    plt.plot(t2[:-1]/1e6, M_acc2[:-1], linestyle="-", linewidth=lw, color="orange", label=r"$\dot{M}_{basic}$, $\alpha$ = 10$^{-3}$")
                    plt.plot(t2[:-1]/1e6, M_acc_max2[:-1], linestyle="--", linewidth=lw, color="orange", label=r"$\dot{M}_{lim}$, $\alpha$ = 10$^{-3}$")
                
                else:
                    plt.plot(t0[:-1]/1e6, M_acc0[:-1], linestyle="-", linewidth=lw, color="black", label="$\dot{M}_{basic}$")
                    plt.plot(t0[:-1]/1e6, M_acc_max0[:-1], linestyle="--", linewidth=lw, color="black", label="$\dot{M}_{lim}$")
                    
                    plt.plot(t0[:-1]/1e6, M_acc0[:-1], linestyle="-", linewidth=lw, color="blue")
                    plt.plot(t0[:-1]/1e6, M_acc_max0[:-1], linestyle="--", linewidth=lw, color="blue")
                    plt.plot(t1[:-1]/1e6, M_acc1[:-1], linestyle="-", linewidth=lw, color="red")
                    plt.plot(t1[:-1]/1e6, M_acc_max1[:-1], linestyle="--", linewidth=lw, color="red")
                    plt.plot(t2[:-1]/1e6, M_acc2[:-1], linestyle="-", linewidth=lw, color="orange")
                    plt.plot(t2[:-1]/1e6, M_acc_max2[:-1], linestyle="--", linewidth=lw, color="orange")
                    
                    plt.legend(loc=(0.1,0.4))
                                
                """
                else:
                    #plt.plot(t0[:-1]/1e6, M_acc_true0[:-1], linestyle="-", linewidth=lw, color="blue", label=r"$\dot{M}$, $\kappa$ = $\kappa_0$")
                    plt.plot(t0[:-1]/1e6, M_acc0[:-1], linestyle="-", linewidth=lw, color="blue", label=r"$\dot{M}_{basic}$, $\kappa$ = $\kappa_0$")
                    plt.plot(t0[:-1]/1e6, M_acc_max0[:-1], linestyle="--", linewidth=lw, color="blue", label=r"$\dot{M}_{lim}$, $\kappa$ = $\kappa_0$")
                    #plt.plot(t1[:-1]/1e6, M_acc_true1[:-1], linestyle="-", linewidth=lw, color="red", label=r"$\dot{M}$, $\kappa$ = $\kappa_0$/6")
                    plt.plot(t1[:-1]/1e6, M_acc1[:-1], linestyle="-", linewidth=lw, color="red", label=r"$\dot{M}_{basic}$, $\kappa$ = $\kappa_0$/6")
                    plt.plot(t1[:-1]/1e6, M_acc_max1[:-1], linestyle="--", linewidth=lw, color="red", label=r"$\dot{M}_{lim}$, $\kappa$ = $\kappa_0$/6")
                    #plt.plot(t2[:-1]/1e6, M_acc_true2[:-1], linestyle="-", linewidth=lw, color="orange", label=r"$\dot{M}$, $\kappa$ = $\kappa_0$/12")
                    plt.plot(t2[:-1]/1e6, M_acc2[:-1], linestyle="-", linewidth=lw, color="orange", label=r"$\dot{M}_{basic}$, $\kappa$ = $\kappa_0$/12")
                    plt.plot(t2[:-1]/1e6, M_acc_max2[:-1], linestyle="--", linewidth=lw, color="orange", label=r"$\dot{M}_{lim}$, $\kappa$ = $\kappa_0$/12")
                """
                
                
                plt.yscale("log")
                plt.ylabel("Accretion Rate / $M_{\oplus}\cdot$yr$^{-1}$")
                #plt.legend()
                #plt.title("Accretion Rate vs. Time")
            #plt.xscale("log")
            #plt.legend()
            ax = plt.gca()
            ax.set_box_aspect(1)
            ax = plt.gca()
            #ax.set_xticks([0,1,2,3])
            #ax.set_xticklabels([0,1,2,3])
            plt.tight_layout()
            plt.show()
        
        
        
        else:
            tau = 1
            
            L, M, t, varis, Mdots = Simulation_Data[0]
            M_acc, M_acc_max, M_acc_true = Mdots
            colors = ["blue","red","orange"]
            
            Figsize, lw, Ticksize, Labelsize, Legendfontsize, Titlesize = Params
            
            plt.rc("figure", figsize=(Figsize,Figsize))
            plt.xlabel("Radius / $R_{J}$")
            plt.rc("xtick", labelsize=Ticksize)
            plt.rc("ytick", labelsize=Ticksize)
            plt.rc("axes", labelsize=Labelsize)
            plt.rc("legend", fontsize=Legendfontsize)
            
            ks = np.trunc(np.linspace(int(t[0]/1e5),int(t[-1]/1e5),len(colors)))
            
            
            import scipy.constants as spc
            G = spc.G * 1e3
            
            lo = 1
            hi = -2
            
            middle_time = 250
            
            m_plot1 = varis[0][0][lo:hi]
            m_plot2 = varis[middle_time][0][lo:hi]
            m_plot3 = varis[-1][0][lo:hi]
            r_plot1 = varis[0][1][lo:hi]
            r_plot2 = varis[middle_time][1][lo:hi]
            r_plot3 = varis[-1][1][lo:hi]
            p_plot1 = varis[0][2][lo:hi]
            p_plot2 = varis[middle_time][2][lo:hi]
            p_plot3 = varis[-1][2][lo:hi]
            k_plot1 = varis[0][3][lo:hi]
            k_plot2 = varis[middle_time][3][lo:hi]
            k_plot3 = varis[-1][3][lo:hi]
            T_plot1 = varis[0][4][lo:hi]
            T_plot2 = varis[middle_time][4][lo:hi]
            T_plot3 = varis[-1][4][lo:hi]
            
            tau_i1 = TLF.Tau_Line_Finder(r_plot1, p_plot1, k_plot1, tau)
            tau_i2 = TLF.Tau_Line_Finder(r_plot2, p_plot2, k_plot2, tau)
            tau_i3 = TLF.Tau_Line_Finder(r_plot3, p_plot3, k_plot3, tau)
            
            
            if y_variable == "T_t,p_t":
                fig, ax1 = plt.subplots()
                color1 = 'red'
                color2 = 'orange'
                
                ax1.set_xlabel('Radius / $R_{J}$')
                ax1.set_ylabel('Density / g$\cdot$cm$^{-3}$', color=color1)
                ax1.loglog(r_plot1/R_J, p_plot1, color=color1, linewidth=lw)
                ax1.loglog(r_plot3/R_J, p_plot3, color=color1, linewidth=lw, linestyle="-.")
                ax1.tick_params(axis='y')
                
                ax2 = ax1.twinx()
                ax2.set_ylabel('Temperature / K', color=color2)
                l1 = ax2.loglog(r_plot1/R_J, T_plot1, color=color2, linewidth=lw, label="0.0 Myr")
                l2 = ax2.loglog(r_plot3/R_J, T_plot3, color=color2, linewidth=lw, linestyle="-.", label="2.0 Myr")
                ax2.tick_params(axis='y')
                
                leg = l1 + l2
                labs = [l.get_label() for l in leg]
                plt.legend(leg, labs, loc=0)
                fig.tight_layout()
            
            
            if y_variable == "p_t comp":
                plt.loglog(r_plot1/R_J, p_plot1, linewidth=lw, linestyle="-", color="blue", label=f"{t[0]/1e6} Myr")
                plt.axvline(x=r_plot1[-tau_i1]/R_J, linewidth=lw, linestyle="--", color="blue")
                plt.loglog(r_plot2/R_J, p_plot2, linewidth=lw, linestyle="-", color="red", label=f"{t[middle_time]/1e6} Myr")
                plt.axvline(x=r_plot2[-tau_i2]/R_J, linewidth=lw, linestyle="--", color="red")
                plt.loglog(r_plot3/R_J, p_plot3, linewidth=lw, linestyle="-", color="orange", label=f"{t[-1]/1e6} Myr")
                plt.axvline(x=r_plot3[-tau_i3]/R_J, linewidth=lw, linestyle="--", color="orange")
                plt.ylabel("Density / g$\cdot$cm$^{-3}$")
                plt.legend(loc="lower left")
            
            if y_variable == "m_t comp":
                m_found1 = 4/3*np.pi*r_plot1**3 * p_plot1 / 500
                m_found2 = 4/3*np.pi*r_plot2**3 * p_plot2 / 50
                m_found3 = 4/3*np.pi*r_plot3**3 * p_plot3 / 50
                #plt.plot(r_plot1/R_J, m_plot1/M_E, linewidth=lw, linestyle="-", color="blue", label=f"{t[0]/1e6} Myr")
                plt.plot(r_plot1/R_J, m_found1/M_E, linewidth=lw, linestyle="-", color="blue", label=f"{t[0]/1e6} Myr")
                #plt.axvline(x=r_plot1[-tau_i1]/R_J, linewidth=lw, linestyle="--", color="blue")
                #plt.plot(r_plot2/R_J, 11*m_plot2/M_E, linewidth=lw, linestyle="-", color="red", label=f"{t[middle_time]/1e6} Myr")
                plt.plot(r_plot2/R_J, m_found2/M_E, linewidth=lw, linestyle="-", color="red", label=f"{t[middle_time]/1e6} Myr")
                #plt.axvline(x=r_plot2[-tau_i2]/R_J, linewidth=lw, linestyle="--", color="red")
                #plt.plot(r_plot3/R_J, 11*m_plot3/M_E, linewidth=lw, linestyle="-", color="orange", label=f"{t[-1]/1e6} Myr")
                plt.plot(r_plot3/R_J, m_found3/M_E, linewidth=lw, linestyle="-", color="orange", label=f"{t[-1]/1e6} Myr")
                #plt.axvline(x=r_plot3[-tau_i3]/R_J, linewidth=lw, linestyle="--", color="orange")
                plt.xscale("log")
                plt.ylabel("Shell Mass / $M_\oplus$")
                plt.legend()
            
            if y_variable == "T_t":
                plt.loglog(r_plot1/R_J, T_plot1, linewidth=lw, linestyle="-", color="blue", label=f"{t[0]/1e6} Myr")
                plt.axvline(x=r_plot1[-tau_i1]/R_J, linewidth=lw, linestyle="--", color="blue")
                plt.loglog(r_plot1[-tau_i1]/R_J,T_plot1[-tau_i1],"x",color="blue",ms=10,mew=3)
                plt.loglog(r_plot2/R_J, T_plot2, linewidth=lw, linestyle="-", color="red", label=f"{t[middle_time]/1e6} Myr")
                plt.axvline(x=r_plot2[-tau_i2]/R_J, linewidth=lw, linestyle="--", color="red")
                plt.loglog(r_plot2[-tau_i2]/R_J,T_plot2[-tau_i2],"x",color="red",ms=10,mew=3)
                plt.loglog(r_plot3/R_J, T_plot3, linewidth=lw, linestyle="-", color="orange", label=f"{t[-1]/1e6} Myr")
                plt.axvline(x=r_plot3[-tau_i3]/R_J, linewidth=lw, linestyle="--", color="orange")
                plt.loglog(r_plot3[-tau_i3]/R_J,T_plot3[-tau_i3],"x",color="orange",ms=10,mew=3)
                plt.ylabel("Temperature / K")
                plt.legend(loc="lower left")
            
            if y_variable == "T_t,p_t diff":
                fig, ax1 = plt.subplots()
                color1 = 'red'
                color2 = 'orange'
                
                ax1.set_xlabel('Radius / $R_{J}$')
                ax1.set_ylabel('Density / g$\cdot$cm$^{-3}$', color=color1)
                ax1.loglog(r_plot1/R_J, p_plot3-p_plot1, color=color1, linewidth=lw)
                ax1.tick_params(axis='y')
                
                ax2 = ax1.twinx()
                ax2.set_ylabel('Temperature / K', color=color2)
                l1 = ax2.loglog(r_plot1/R_J, T_plot3-T_plot1, color=color2, linewidth=lw)
                ax2.tick_params(axis='y')
                
                leg = l1 + l2
                labs = [l.get_label() for l in leg]
                plt.legend(leg, labs, loc=0)
                fig.tight_layout()
            
            
            for j in range(len(ks)):
                
                tk = int(ks[j]*1e5)
                k = np.where(t == tk)[0][0]
                
                #lo = varis[k][1]
                #hi = varis[k][2]
                
                m_plot = varis[k][0][lo:hi]
                r_plot = varis[k][1][lo:hi]
                p_plot = varis[k][2][lo:hi]
                k_plot = varis[k][3][lo:hi]
                T_plot = varis[k][4][lo:hi]
                
                tau_i = TLF.Tau_Line_Finder(r_plot, p_plot, k_plot, tau)
                phi_plot = -(G*(m_plot+(10*M_E))/r_plot)
                
                if y_variable == "T_t old":
                    plt.loglog(r_plot/R_J, T_plot, linewidth=lw, linestyle="-", color=colors[j], label=f"{t[k]/1e6} Myr")
                    plt.axvline(x=r_plot[-tau_i]/R_J, linewidth=lw, linestyle="--", color=colors[j])
                    plt.loglog(r_plot[-tau_i]/R_J,T_plot[-tau_i],"x",color=colors[j],ms=10,mew=3)
                    plt.ylabel("Temperature / K")
                    #plt.title("Temperature vs. Radius")
                
                elif y_variable == "phi_t":
                    plt.loglog(r_plot/R_J, np.abs(phi_plot), linewidth=lw, linestyle="-", color=colors[j], label=f"{t[k]/1e6} Myr")
                    plt.xscale("log")
                    plt.ylabel("Phi")
                    #plt.title("Phi vs. Radius")
            
            if y_variable == "T_t,p_t":
                ax = plt.gca()
                leg = ax.get_legend()
                leg.legendHandles[0].set_color('black')
                leg.legendHandles[1].set_color('black')
            plt.axvline(r_plot[0]/R_J, linestyle="--", color="black")
            plt.axvline(r_plot[-1]/R_J, linestyle="--", color="black")
            ax = plt.gca()
            ax.set_box_aspect(1)
            ax = plt.gca()
            ax.set_xticks([r_plot[0]/R_J,1,10,r_plot[-1]/R_J])
            ax.set_xticklabels(["$r_{core}$",1,10,"$r_{p}$"])
            plt.tight_layout()
            plt.show()
    
    
    
    else:
        Simulation_Data = Simulation_Data[0]
        var_plot = Simulation_Data._var.copy()
        if Movshovitz == True:
            m, r, p, k, T, P, s = [var[-Simulation_Data._break:Simulation_Data._hi] for var in var_plot]
        else:
            m, r, p, k, T, P, s = [var[Simulation_Data._lo:Simulation_Data._hi] for var in var_plot]
        
        for x in range(len(Simulation_Data._boundary_i)):
            print(r[Simulation_Data._boundary_i[x]])
        
        if type(Params) != type(None):
            Figsize, lw, Ticksize, Labelsize, Legendfontsize, Titlesize = Params
        
            plt.rc("figure", figsize=(Figsize,Figsize))
            plt.xlabel("Radius / $R_{J}$")
            plt.rc("xtick", labelsize=Ticksize)
            plt.rc("ytick", labelsize=Ticksize)
            plt.rc("axes", labelsize=Labelsize)
            plt.rc("legend", fontsize=Legendfontsize)
        
        # Plots against m
        if x_variable == "m":
            plt.xlabel("Mass / $M_{\oplus}$")
            
            if y_variable == "r":
                plt.loglog(m/M_E, r/R_J, linewidth=lw)
                plt.ylabel("Radius / $R_{J}$")
                #plt.title("Radius vs. Mass")
            
            elif y_variable == "p":
                plt.plot(m/M_E, p, linewidth=lw)
                plt.ylabel("Density / g$\cdot$cm$^{-3}$")
                #plt.title("Density vs. Mass")
            
            elif y_variable == "k":
                plt.plot(m/M_E, k, linewidth=lw)
                plt.ylabel("Opacity / cm$^{2}$$\cdot$g$^{-1}$")
                #plt.title("Opacity vs. Mass")
            
            elif y_variable == "T":
                plt.plot(m/M_E, T, linewidth=lw)
                plt.ylabel("Temperature / K")
                #plt.title("Temperature vs. Mass")
            
            elif y_variable == "P":
                plt.plot(m/M_E, P, linewidth=lw)
                plt.ylabel("Pressure / dyne$\cdot$cm$^{-2}$")
                #plt.title("Pressure vs. Mass")
            
            elif y_variable == "s":
                plt.plot(m/M_E, s, linewidth=lw)
                plt.ylabel(r"Measure of Entropy ($T/P^{(\gamma-1)/\gamma}$)")
                #plt.title("Measure of Entropy vs. Mass")
        
        # Plots against r
        elif x_variable == "r":
            if y_variable != "p+T" and y_variable != "k+p+T":
                plt.xlabel("Radius / $R_{J}$")
                
                if y_variable == "m":
                    plt.loglog(r/R_J, m/M_E, linewidth=lw)
                    #plt.loglog(r/R_J, m/M_E, "x")
                    plt.ylabel("Mass / $M_{\oplus}$")
                    #plt.title("Mass vs. Radius")
                
                elif y_variable == "p":
                    plt.loglog(r/R_J, p, linewidth=lw, color="blue", label="Our Simulation")
                    #plt.loglog(r/R_J, p, "x", label="Our Simulation")
                    if Comparison_Data != None:
                        plt.loglog(Comparison_Data[0][0]/R_J, Comparison_Data[0][1], "--", linewidth=lw, color="orange", label="M10")
                        plt.legend()
                    plt.ylabel("Density / g$\cdot$cm$^{-3}$")
                    #plt.title("Density vs. Radius")
                
                elif y_variable == "k":
                    plt.loglog(r/R_J, k, linewidth=lw, color="blue", label="Our Simulation")
                    #plt.loglog(r/R_J, k, "x", label="Our Simulation")
                    if Comparison_Data != None:
                        plt.loglog(Comparison_Data[1][0]/R_J, Comparison_Data[1][1], "--", linewidth=lw, color="orange", label="M10")
                        plt.legend()
                    plt.ylabel("Opacity / cm$^{2}$$\cdot$g$^{-1}$")
                    #plt.title("Opacity vs. Radius")
    
                
                elif y_variable == "T":
                    plt.loglog(r/R_J, T, linewidth=lw, color="blue", label="Our Simulation")
                    #plt.plot(r/R_J, T, "x", label="Our Simulation")
                    if Comparison_Data != None:
                        plt.loglog(Comparison_Data[2][0]/R_J, Comparison_Data[2][1], "--", linewidth=lw, color="orange", label="M10")
                        plt.legend()
                    plt.ylabel("Temperature / K")
                    #plt.title("Temperature vs. Radius")
                
                elif y_variable == "P":
                    plt.loglog(r/R_J, P, linewidth=lw)
                    #plt.plot(r/R_J, P, "x")
                    plt.ylabel("Pressure / dyne$\cdot$cm$^{-2}$")
                    #plt.title("Pressure vs. Radius")
                
                elif y_variable == "s":
                    plt.loglog(r/R_J, s, linewidth=lw, color="blue")
                    #plt.loglog(r/R_J, s, "x")
                    plt.ylabel(r"Measure of Entropy ($T/P^{(\gamma-1)/\gamma}$)")
                    #plt.title("Measure of Entropy vs. Radius")
            
            elif y_variable == "p+T":
                fig, ax1 = plt.subplots()
                color1 = 'red'
                color2 = 'orange'
                color3 = 'blue'
                color4 = 'green'
                
                ax1.set_xlabel('Radius / $R_{J}$')
                ax1.set_ylabel('Density / g$\cdot$cm$^{-3}$', color=color1)
                l1 = ax1.loglog(r/R_J, p, color=color1, linewidth=lw)
                if Comparison_Data != None:
                    l3 = ax1.loglog(Comparison_Data[0][0]/R_J, Comparison_Data[0][1], color=color3, linestyle="--", label="Movshovitz Density")
                ax1.tick_params(axis='y')
                
                ax2 = ax1.twinx()
                ax2.set_ylabel('Temperature / K', color=color2)
                l2 = ax2.loglog(r/R_J, T, color=color2, linewidth=lw)
                if Comparison_Data != None:
                    l4 = ax2.loglog(Comparison_Data[2][0]/R_J, Comparison_Data[2][1], color=color4, linestyle="--", label="Movshovitz Temperature")
                ax2.tick_params(axis='y')
                fig.tight_layout()
                
                if Comparison_Data != None:
                    leg = l1 + l3 + l2 + l4
                else:
                    leg = l1 + l2
                labs = [l.get_label() for l in leg]
                #plt.legend(leg, labs, loc=0)
                
            
            elif y_variable == "k+p+T":
                
                for i in range(1,len(T)+1):
                    boundary_k1 = i
                    if T[-i] > 4000:
                        break
                
                for i in range(1,len(T)+1):
                    boundary_k2 = i
                    if P[i] < 1:
                        break
                
                fig, ax1 = plt.subplots()
                color1 = 'blue'
                color2 = 'red'
                color3 = 'orange'
                
                ax1.set_xlabel('Radius / $R_{J}$')
                ax1.set_ylabel('Density / g$\cdot$cm$^{-3}$', color=color2)
                ax1.loglog(r/R_J, p, color=color2, linewidth=lw)
                ax1.tick_params(axis='y')
                
                ax2 = ax1.twinx()
                ax2.set_ylabel('Temperature / K', color=color3)
                ax2.loglog(r/R_J, T, color=color3, linewidth=lw)
                ax2.tick_params(axis='y')
                fig.tight_layout()
                
                ax3 = ax1.twinx()
                ax3.set_ylabel('Opacity / cm$^{2}$$\cdot$g$^{-1}$', color=color1)
                ax3.spines.right.set_position(("axes", 1.3))
                ax3.loglog(r[:-boundary_k1]/R_J, k[:-boundary_k1], color=color1, linewidth=lw, linestyle="--")
                ax3.loglog(r[-boundary_k1:boundary_k2]/R_J, k[-boundary_k1:boundary_k2], color=color1, linewidth=lw)
                ax3.loglog(r[boundary_k2:]/R_J, k[boundary_k2:], color=color1, linewidth=lw, linestyle="--")
                ax3.tick_params(axis='y')
                fig.tight_layout()
            
            
            elif y_variable == "dT4_dm":
                dT4_dm = []
                r_plot = []
                for i in range(len(r)-1):
                    if np.isfinite(r[i]):
                        dT4_dm.append(abs(((T[i]**4)-(T[i+1]**4))/(m[i]-m[i+1])))
                        r_plot.append(r[i])
                plt.loglog(r_plot/R_J, dT4_dm)
                plt.ylabel("abs(dT4_dm)")
                #plt.title("dT4_dm vs. Radius")
        
        #for x in range(len(Simulation_Data._boundary_i)):
        #    plt.axvline(x=r[-Simulation_Data._boundary_i[x][0]]/R_J, linestyle="--", color="black")
        darkness = 0.2
        colour = "black"
        if Comparison_Data != None:
            rcore = Comparison_Data[0][0][0]/R_J
        else:
            rcore = r[0]/R_J
        if len(Simulation_Data._boundary_i) != 0:
            if Movshovitz==True:
                plt.axvspan(r[0]/R_J, r[-Simulation_Data._boundary_i[3][0]]/R_J, alpha=darkness, facecolor=colour)
                plt.axvspan(r[-Simulation_Data._boundary_i[2][0]]/R_J, r[-Simulation_Data._boundary_i[1][0]]/R_J, alpha=darkness, facecolor=colour)
                plt.axvspan(r[-Simulation_Data._boundary_i[0][0]]/R_J, r[-1]/R_J, alpha=darkness, facecolor=colour)
            else:
                plt.axvspan(r[0]/R_J, r[-Simulation_Data._boundary_i[0][0]]/R_J, alpha=darkness, facecolor=colour)
        plt.axvline(rcore, linestyle="--", color="black")
        plt.axvline(r[-1]/R_J, linestyle="--", color="black")
        ax = plt.gca()
        ax.set_box_aspect(1)
        ax = plt.gca()
        ax.set_xticks([rcore,1,10,r[-1]/R_J])
        ax.set_xticklabels(["$r_{core}$",1,10,"$r_{p}$"])
        plt.tight_layout()
        plt.show()
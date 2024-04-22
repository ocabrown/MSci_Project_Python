#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 23:04:24 2022

@author: ollie
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import Mass_Grid_Solver as MGS
import IsoT_Model as ITM
import IsoD_Model as IDM
import IsoO_Model as IOM
import GPF_Movshovitz_Opacity as GPFMO
import GPF_Freedman_Opacity as GPFFO
import GPF_Combined_Opacity as GPFCO
import Maximum_Accretion_Rate_Calculator as MARC
import L_Finder_IsoO as LFIO
import L_Finder_FreO as LFFO
import L_Finder_ComO as LFCO


def N1(x,a):
    return -1*x + a

def N2(x,a):
    return -2*x + a


"""
from decimal import Decimal
F = 1e16
kappa = 1e-18
alpha = 1e-13
ntot_init = 1e3

lhs = Decimal(F)*Decimal(kappa)
a = Decimal(alpha)/lhs
b = Decimal(1)
c = -Decimal(ntot_init)
nplus = (-b + np.sqrt(b*b-Decimal(4)*a*c))/Decimal(2.)/a
fplus = nplus / Decimal(ntot_init)
fzero = (Decimal(1)-fplus)


F = 1e16
kappa = 1e-18
alpha = 1e-13
ntot_init = 1e3

lhs = F*kappa
a = alpha/lhs
b = 1
c = -ntot_init
nplus = (-b + np.sqrt(b*b-4*a*c))/2./a
fplus = nplus / ntot_init
fzero = (1-fplus)
"""



# Isothermal Convergence Function

def ConvergenceIsoT(mc, me, rp, pd, Td, n, acc_m, p_ana_c_ini, Plot=True):
    
    def straight_line(x, m, c):
        return m * x + c
    
    L2norms_p = []
    L2norms_dP_dr = []
    min_dP_drs = []
    
    for i in range(len(n)):
        m_grid_isoT, _ = MGS.Mass_Grid_Solver_Iso_Temperature(mc, me, rp, pd, Td, n[i], acc_m, p_ana_c_ini)
        IsoT = ITM.Iso_Temperature_Giant_Planet(mc,me,rp,pd,Td,n[i],grid_given=True,m_grid=m_grid_isoT)
        IsoT.Planet_Formation()
        
        L2norms_p.append(IsoT.L2norm()[0])
        L2norms_dP_dr.append(IsoT.L2norm()[1])
        min_dP_drs.append(IsoT.L2norm()[2])
    
    L2norms_p = np.array(L2norms_p)
    L2norms_dP_dr = np.array(L2norms_dP_dr)
    min_dP_drs = np.array(min_dP_drs)
    
    vals1, val_errs1 = curve_fit(straight_line, np.log10(n), np.log10(L2norms_p))
    order1, order_std1 = vals1[0], np.sqrt(val_errs1[0][0])
    print("Order =", order1, "±", order_std1)
    
    vals2, val_errs2 = curve_fit(straight_line, np.log10(n), np.log10(L2norms_dP_dr))
    order2, order_std2 = vals2[0], np.sqrt(val_errs2[0][0])
    print("Order =", order2, "±", order_std2)
    
    p_popt1, p_pcov1 = curve_fit(N1,np.log10(n), np.log10(L2norms_p))
    P_popt1, P_pcov1 = curve_fit(N1,np.log10(n), np.log10(L2norms_dP_dr))
    p_popt2, p_pcov2 = curve_fit(N2,np.log10(n), np.log10(L2norms_p))
    P_popt2, P_pcov2 = curve_fit(N2,np.log10(n), np.log10(L2norms_dP_dr))
    
    if n[-1] > 1e3:
        ind = np.where(n==1e3)
        p_popt1, p_pcov1 = curve_fit(N1,np.log10(n[:int(ind[0][0]-5)]), np.log10(L2norms_p[:int(ind[0][0]-5)]))
        p_popt2, p_pcov2 = curve_fit(N2,np.log10(n[:int(ind[0][0]-5)]), np.log10(L2norms_p[:int(ind[0][0]-5)]))
    
    
    Figsize, lw, Ticksize, Labelsize, Legendfontsize = [7, 2, 20, 20, 15]
    
    plt.rc("figure", figsize=(Figsize,Figsize))
    plt.rc("xtick", labelsize=Ticksize)
    plt.rc("ytick", labelsize=Ticksize)
    plt.rc("axes", labelsize=Labelsize)
    plt.rc("legend", fontsize=Legendfontsize)
    
    if Plot == True:
        plt.figure()
        plt.loglog(n,L2norms_p,color="blue",linewidth=lw)
        plt.loglog(n,10**N1(np.log10(n),p_popt1),color="orange",linestyle="--",linewidth=lw)
        plt.loglog(n,10**N2(np.log10(n),p_popt2),color="red",linestyle="--",linewidth=lw)
        plt.xlabel("Number of Cells")
        plt.ylabel("L$^2$-norm")
        #plt.title("Density L$^2$ Norm vs. Number of Cells")
        ax = plt.gca()
        ax.set_box_aspect(1)
        plt.tight_layout()
        plt.show()
        
        x=1.6
        
        plt.figure()
        plt.loglog(n,L2norms_dP_dr*x,color="blue",linewidth=lw)
        plt.loglog(n,10**N1(np.log10(n),P_popt1)*x,color="orange",linestyle="--",linewidth=lw)
        plt.loglog(n,10**N2(np.log10(n),P_popt2)*x,color="red",linestyle="--",linewidth=lw)
        plt.xlabel("Number of Cells")
        plt.ylabel("L$^2$-norm")
        #plt.title("Pressure Gradient L$^2$ Norm vs. Number of Cells")
        ax = plt.gca()
        ax.set_box_aspect(1)
        plt.tight_layout()
        plt.show()
        
        plt.figure()
        plt.loglog(n, abs(min_dP_drs),linewidth=lw)
        plt.xlabel("Number of Cells")
        plt.ylabel("Minimum Pressure Gradient")
        #plt.title("Minimum Pressure Gradient vs. Number of Cells")
        ax = plt.gca()
        ax.set_box_aspect(1)
        plt.tight_layout()
        plt.show()
        
        plt.figure()
        plt.loglog(n,L2norms_p,linewidth=lw,color="red",label="Density L$^2$-norm")
        plt.loglog(n,L2norms_dP_dr,linewidth=lw,color="blue",label="dP/dr L$^2$-norm")
        plt.loglog(n,10**N1(np.log10(n),P_popt1),color="orange",linestyle="--",linewidth=lw)
        plt.xlabel("Number of Cells")
        plt.ylabel("L$^2$-norm")
        ax = plt.gca()
        ax.set_box_aspect(1)
        plt.tight_layout()
        plt.legend()
        plt.show()
        
        """
        fig, ax1 = plt.subplots()
        color1 = 'tab:red'
        ax1.set_xlabel('Number of Cells')
        ax1.set_ylabel("Density L$^2$ Norm", color=color1)
        l1 = ax1.loglog(n, L2norms_p, color=color1, linewidth=lw, label="Density L$^2$ Norm")
        ax1.tick_params(axis='y', labelcolor=color1)
        
        ax2 = ax1.twinx()
        
        color1 = 'tab:orange'
        ax2.set_ylabel('Pressure Gradient L$^2$ Norm', color=color1)
        l2 = ax2.loglog(n, L2norms_dP_dr, color=color1, linewidth=lw, label="Pressure Gradient L$^2$ Norm")
        ax2.tick_params(axis='y', labelcolor=color1)
        
        ax = plt.gca()
        ax.set_box_aspect(1)
        plt.tight_layout()
        
        leg = l1 + l2
        labs = [l.get_label() for l in leg]
        plt.legend(leg, labs, loc=9)
        """
    
    return L2norms_p, L2norms_dP_dr



# Constant Density Convergence Function

def ConvergenceIsoD(mc, me, rp, pd, kd, Td, lp, n, Plot=True):
    
    def straight_line(x, m, c):
        return m * x + c
    
    L2norms_T = []
    L2norms_dP_dr = []
    L2norms_dP_dr_lower = []
    L2norms_dP_dr_upper = []
    min_dP_drs = []
    
    for i in range(len(n)):
        IsoD = IDM.Iso_Density_Giant_Planet(mc,me,rp,pd,kd,Td,lp,n[i])
        IsoD.Planet_Formation()
        
        L2norms_T.append(IsoD.L2norm()[0])
        L2norms_dP_dr.append(IsoD.L2norm()[1])
        L2norms_dP_dr_lower.append(IsoD.L2norm()[2])
        L2norms_dP_dr_upper.append(IsoD.L2norm()[3])
        min_dP_drs.append(IsoD.L2norm()[4])
    
    L2norms_T = np.array(L2norms_T)
    L2norms_dP_dr = np.array(L2norms_dP_dr)
    L2norms_dP_dr_lower = np.array(L2norms_dP_dr_lower)
    L2norms_dP_dr_upper = np.array(L2norms_dP_dr_upper)
    min_dP_drs = np.array(min_dP_drs)
    
    vals, val_errs = curve_fit(straight_line, np.log10(n), np.log10(L2norms_T))
    order, order_std = vals[0], np.sqrt(val_errs[0][0])
    print(order, "±", order_std)
    
    T_popt1, T_pcov1 = curve_fit(N1,np.log10(n), np.log10(L2norms_T))
    P_popt1, P_pcov1 = curve_fit(N1,np.log10(n), np.log10(L2norms_dP_dr_upper))
    T_popt2, T_pcov2 = curve_fit(N2,np.log10(n), np.log10(L2norms_T))
    P_popt2, P_pcov2 = curve_fit(N2,np.log10(n), np.log10(L2norms_dP_dr_upper))
    
    Figsize, lw, Ticksize, Labelsize, Legendfontsize = [7, 2, 20, 20, 15]
    
    plt.rc("figure", figsize=(Figsize,Figsize))
    plt.rc("xtick", labelsize=Ticksize)
    plt.rc("ytick", labelsize=Ticksize)
    plt.rc("axes", labelsize=Labelsize)
    plt.rc("legend", fontsize=Legendfontsize)
    
    if Plot == True:
        plt.figure()
        plt.loglog(n,L2norms_T, linewidth=lw)
        plt.loglog(n,10**N1(np.log10(n),T_popt1),color="orange",linestyle="--",linewidth=lw)
        plt.loglog(n,10**N2(np.log10(n),T_popt2),color="red",linestyle="--",linewidth=lw)
        plt.xlabel("Number of Cells")
        plt.ylabel("Temperature L$^2$ Norm")
        #plt.title("Temperature L$^2$ Norm vs. Number of Cells")
        ax = plt.gca()
        ax.set_box_aspect(1)
        plt.tight_layout()
        plt.show()
        
        plt.figure()
        plt.plot(n,L2norms_dP_dr, linewidth=lw)
        plt.xlabel("Number of Cells")
        plt.ylabel("Pressure Gradient L$^2$ Norm")
        #plt.title("Pressure Gradient L$^2$ Norm vs. Number of Cells")
        ax = plt.gca()
        ax.set_box_aspect(1)
        plt.tight_layout()
        plt.show()
        
        plt.figure()
        plt.plot(n,L2norms_dP_dr_lower, linewidth=lw)
        plt.xlabel("Number of Cells")
        plt.ylabel("Pressure Gradient L$^2$ Norm")
        #plt.title("Lower Pressure Gradient L$^2$ Norm vs. Number of Cells")
        ax = plt.gca()
        ax.set_box_aspect(1)
        plt.tight_layout()
        plt.show()
        
        plt.figure()
        plt.loglog(n,L2norms_dP_dr_upper, linewidth=lw)
        plt.loglog(n,10**N1(np.log10(n),P_popt1),color="orange",linestyle="--",linewidth=lw)
        plt.loglog(n,10**N2(np.log10(n),P_popt2),color="red",linestyle="--",linewidth=lw)
        plt.xlabel("Number of Cells")
        plt.ylabel("Pressure Gradient L$^2$ Norm")
        #plt.title("Upper Pressure Gradient L$^2$ Norm vs. Number of Cells")
        ax = plt.gca()
        ax.set_box_aspect(1)
        plt.tight_layout()
        plt.show()
        
        plt.figure()
        plt.plot(n, abs(min_dP_drs), linewidth=lw)
        plt.xlabel("Number of Cells")
        plt.ylabel("Minimum Pressure Gradient")
        #plt.title("Minimum Pressure Gradient vs. Number of Cells")
        ax = plt.gca()
        ax.set_box_aspect(1)
        plt.tight_layout()
        plt.show()
        
        plt.figure()
        plt.loglog(n,L2norms_T,linewidth=lw,color="red",label="Temperature L$^2$-norm")
        plt.loglog(n,10**N1(np.log10(n),T_popt1),color="orange",linestyle="--",linewidth=lw)
        plt.loglog(n,L2norms_dP_dr,linewidth=lw,color="blue",label="dP/dr L$^2$-norm")
        plt.loglog(n,10**N1(np.log10(n),P_popt1),color="orange",linestyle="--",linewidth=lw)
        plt.xlabel("Number of Cells")
        plt.ylabel("L$^2$-norm")
        ax = plt.gca()
        ax.set_box_aspect(1)
        plt.tight_layout()
        plt.legend()
        plt.show()
        
        fig, ax1 = plt.subplots()
        color1 = 'tab:red'
        ax1.set_xlabel('Number of Cells')
        ax1.set_ylabel("Temperature L$^2$ Norm", color=color1)
        l1 = ax1.loglog(n, L2norms_T, color=color1, linewidth=lw, label="Temperature L$^2$ Norm")
        ax1.tick_params(axis='y', labelcolor=color1)
        
        ax2 = ax1.twinx()
        
        color1 = 'tab:orange'
        ax2.set_ylabel('Pressure Gradient L$^2$ Norm', color=color1)
        l2 = ax2.loglog(n, L2norms_dP_dr_upper, color=color1, linewidth=lw, label="Pressure Gradient L$^2$ Norm")
        ax2.tick_params(axis='y', labelcolor=color1)
        
        ax = plt.gca()
        ax.set_box_aspect(1)
        plt.tight_layout()
        
        leg = l1 + l2
        labs = [l.get_label() for l in leg]
        plt.legend(leg, labs, loc=9)
    
    return L2norms_T, L2norms_dP_dr_upper



# Constant Opacity Convergence Function

def ConvergenceIsoO(mc, me, rp, pd, kd, Td, lp, n, acc_m, p_ana_c_ini, Plot=True):
    
    def straight_line(x, m, c):
        return m * x + c
    
    r0s = []
    L2norms_dP_dr = []
    min_dP_drs = []
    
    m_grid_isoO, _ = MGS.Mass_Grid_Solver_Iso_Opacity(mc, me, rp, pd, kd, Td, lp, 10*n[-1], acc_m, p_ana_c_ini)
    IsoO = IOM.Iso_Opacity_Giant_Planet(mc, me, rp, pd, kd, Td, lp, 10*n[-1], grid_given=True, m_grid=m_grid_isoO)
    r0_high = IsoO.Planet_Formation()[1][1]
    
    for i in range(len(n)):
        m_grid_isoO, _ = MGS.Mass_Grid_Solver_Iso_Opacity(mc, me, rp, pd, kd, Td, lp, n[i], acc_m, p_ana_c_ini)
        IsoO = IOM.Iso_Opacity_Giant_Planet(mc, me, rp, pd, kd, Td, lp, n[i], grid_given=True, m_grid=m_grid_isoO)
        IsoO_Sim_Data = IsoO.Planet_Formation()
        r0s.append(IsoO_Sim_Data[1][1])
        
        L2norms_dP_dr.append(IsoO.L2norm()[0])
        min_dP_drs.append(IsoO.L2norm()[1])
    
    L2norms_dP_dr = np.array(L2norms_dP_dr)
    min_dP_drs = np.array(min_dP_drs)
    
    vals2, val_errs2 = curve_fit(straight_line, np.log10(n), np.log10(L2norms_dP_dr))
    order2, order_std2 = vals2[0], np.sqrt(val_errs2[0][0])
    print("Order =", order2, "±", order_std2)
    
    P_popt1, P_pcov1 = curve_fit(N1,np.log10(n), np.log10(L2norms_dP_dr))
    P_popt2, P_pcov2 = curve_fit(N2,np.log10(n), np.log10(L2norms_dP_dr))
    
    r_popt1, r_pcov1 = curve_fit(N1,np.log10(n), np.log10(abs((r0_high-r0s)/r0_high)))
    
    Figsize, lw, Ticksize, Labelsize = [7, 2, 20, 20]
    
    plt.rc("figure", figsize=(Figsize,Figsize))
    plt.rc("xtick", labelsize=Ticksize)
    plt.rc("ytick", labelsize=Ticksize)
    plt.rc("axes", labelsize=Labelsize)
    
    # Convergence Plots
    plt.figure()
    plt.loglog(n,L2norms_dP_dr,linewidth=lw)
    plt.loglog(n,10**N1(np.log10(n),P_popt1),color="orange",linestyle=":",linewidth=lw)
    plt.loglog(n,10**N2(np.log10(n),P_popt2),color="red",linestyle=":",linewidth=lw)
    plt.xlabel("Number of Cells")
    plt.ylabel("Pressure Gradient L$^2$ Norm")
    #plt.title("Pressure Gradient L$^2$ Norm vs. Number of Cells")
    ax = plt.gca()
    ax.set_box_aspect(1)
    plt.tight_layout()
    plt.show()
    
    plt.figure()
    plt.loglog(n, abs(min_dP_drs),linewidth=lw)
    plt.xlabel("Number of Cells")
    plt.ylabel("Minimum Pressure Gradient")
    #plt.title("Minimum Pressure Gradient vs. Number of Cells")
    ax = plt.gca()
    ax.set_box_aspect(1)
    plt.tight_layout()
    plt.show()
    
    plt.figure()
    plt.loglog(n, abs((r0_high-r0s)/r0_high),color="blue", linewidth=lw)
    plt.loglog(n,10**N1(np.log10(n),r_popt1),color="orange",linestyle="--",linewidth=lw)
    plt.xlabel("Number of Cells")
    plt.ylabel("$r_{inner}$ Relative Difference")
    #plt.title("$r_{inner}$ Relative Difference vs. Number of Cells")
    ax = plt.gca()
    ax.set_box_aspect(1)
    plt.tight_layout()
    plt.show()
    
    return abs((r0s[-1]-r0s[:-1])/r0s[-1])


# Movshovitz Opacity Convergence Function

def ConvergenceMovO(mc, me, rp, pd, kd, Td, lp, n, Data, Plot=True):
    
    r0s = []
    
    for i in range(len(n)):
        MovO = GPFMO.Movshovitz_Opacity_Giant_Planet(mc, me, rp, pd, kd, Td, lp, n[i])
        MovO_Sim_Data = MovO.Planet_Formation(Data)
        r0s.append(MovO_Sim_Data[1][-(int(MovO._boundary_i[0][0]))])
        #r0s.append(MovO_Sim_Data[1][-(int(MovO._break))])
        #r0s.append(MovO_Sim_Data[1][1])
    
    Figsize, lw, Ticksize, Labelsize = [7, 2, 20, 20]
    
    plt.rc("figure", figsize=(Figsize,Figsize))
    plt.rc("xtick", labelsize=Ticksize)
    plt.rc("ytick", labelsize=Ticksize)
    plt.rc("axes", labelsize=Labelsize)
    
    # Convergence Plot
    plt.figure()
    plt.loglog(n[:-1], abs((r0s[-1]-r0s[:-1])/r0s[-1]), linewidth=lw)
    plt.xlabel("Number of Cells")
    plt.ylabel("$r_{inner}$ Relative Difference")
    #plt.title("$r_{inner}$ Relative Difference vs. Number of Cells")
    ax = plt.gca()
    ax.set_box_aspect(1)
    plt.tight_layout()
    plt.show()
    
    return abs((r0s[-1]-r0s[:-1])/r0s[-1])



# Freedman Opacity Convergence Function

def ConvergenceFreO(mc, me, rp, pd, kd, Td, lp, met, acc_k, n, acc_m, p_ana_c_ini, Plot=True):
    
    r0s = []
    
    #m_grid_fre, _ = MGS.Mass_Grid_Solver_Freedman_Opacity(mc, me, rp, pd, kd, Td, lp, met, acc_k, 10*n[-1], acc_m, p_ana_c_ini)
    FreO = GPFFO.Freedman_Opacity_Giant_Planet(mc, me, rp, pd, kd, Td, lp, met, acc_k, 10*n[-1], p_ana_c=p_ana_c_ini)#, grid_given=True, m_grid=m_grid_fre)
    r0_high = FreO.Planet_Formation()[1][1]
    
    for i in range(len(n)):
        #m_grid_fre, _ = MGS.Mass_Grid_Solver_Freedman_Opacity(mc, me, rp, pd, kd, Td, lp, met, acc_k, n[i], acc_m, p_ana_c_ini)
        FreO = GPFFO.Freedman_Opacity_Giant_Planet(mc, me, rp, pd, kd, Td, lp, met, acc_k, n[i], p_ana_c=p_ana_c_ini)#grid_given=True, m_grid=m_grid_fre)
        FreO_Sim_Data = FreO.Planet_Formation()
        r0s.append(FreO_Sim_Data[1][1])
    
    Figsize, lw, Ticksize, Labelsize = [7, 2, 20, 20]
    
    plt.rc("figure", figsize=(Figsize,Figsize))
    plt.rc("xtick", labelsize=Ticksize)
    plt.rc("ytick", labelsize=Ticksize)
    plt.rc("axes", labelsize=Labelsize)
    
    # Convergence Plots
    plt.figure()
    plt.loglog(n, abs((r0_high-r0s)/r0_high), linewidth=lw)
    plt.xlabel("Number of Cells")
    plt.ylabel("$r_{inner}$ Relative Difference")
    #plt.title("$r_{inner}$ Relative Difference vs. Number of Cells")
    ax = plt.gca()
    ax.set_box_aspect(1)
    plt.tight_layout()
    plt.show()
    
    return abs((r0s[-1]-r0s[:-1])/r0s[-1])



# Freedman Opacity Convergence Function

def ConvergenceComO(mc, me, rp, pd, kd, Td, lp, ad, d_to_g, acc_k, n, acc_m, p_ana_c_ini, Plot=True):
    
    r0s = []
    
    m_grid_com, _ = MGS.Mass_Grid_Solver_Combined_Opacity(mc, me, rp, pd, kd, Td, lp, ad, d_to_g, acc_k, 10*n[-1], acc_m, p_ana_c_ini)
    ComO = GPFFO.Combined_Opacity_Giant_Planet(mc, me, rp, pd, kd, Td, lp, ad, d_to_g, acc_k, 10*n[-1], grid_given=True, m_grid=m_grid_com)
    r0_high = ComO.Planet_Formation()[1][1]
    
    for i in range(len(n)):
        m_grid_com, _ = MGS.Mass_Grid_Solver_Combined_Opacity(mc, me, rp, pd, kd, Td, lp, ad, d_to_g, acc_k, n[i], acc_m, p_ana_c_ini)
        ComO = GPFFO.Combined_Opacity_Giant_Planet(mc, me, rp, pd, kd, Td, lp, ad, d_to_g, acc_k, n[i], grid_given=True, m_grid=m_grid_com)
        ComO_Sim_Data = ComO.Planet_Formation()
        r0s.append(ComO_Sim_Data[1][1])
    
    Figsize, lw, Ticksize, Labelsize = [7, 2, 20, 20]
    
    plt.rc("figure", figsize=(Figsize,Figsize))
    plt.rc("xtick", labelsize=Ticksize)
    plt.rc("ytick", labelsize=Ticksize)
    plt.rc("axes", labelsize=Labelsize)
    
    # Convergence Plots
    plt.figure()
    plt.loglog(n, abs((r0_high-r0s)/r0_high), linewidth=lw)
    #plt.loglog(n[:-1], abs((r0s[-1]-r0s[:-1])/r0s[-1]), linewidth=lw)
    plt.xlabel("Number of Cells")
    plt.ylabel("$r_{inner}$ Relative Difference")
    #plt.title("$r_{inner}$ Relative Difference vs. Number of Cells")
    ax = plt.gca()
    ax.set_box_aspect(1)
    plt.tight_layout()
    plt.show()
    
    return abs((r0s[-1]-r0s[:-1])/r0s[-1])



# Sigma Convergence Function

def ConvergenceSigma(Val_outer, Ms, Mp, Tp, Rp, r_in, r_out, alpha, n, Plot=True):
    
    S0s = []
    
    _, Sigma, _ = MARC.Max_Acc_Rate(Val_outer, Ms, Mp, Tp, Rp, r_in, r_out, alpha, 10*n[-1])
    S0_high = Sigma[0]
    
    for i in range(len(n)):
        _, Sigma, _ = MARC.Max_Acc_Rate(Val_outer, Ms, Mp, Tp, Rp, r_in, r_out, alpha, n[i])
        S0s.append(Sigma[0])
    
    S_popt1, S_pcov1 = curve_fit(N1, np.log10(n), np.log10(abs((S0_high-S0s)/S0_high)))
    
    Figsize, lw, Ticksize, Labelsize = [7, 2, 20, 20]
    
    plt.rc("figure", figsize=(Figsize,Figsize))
    plt.rc("xtick", labelsize=Ticksize)
    plt.rc("ytick", labelsize=Ticksize)
    plt.rc("axes", labelsize=Labelsize)
    
    # Convergence Plots
    plt.figure()
    plt.loglog(n, abs((S0_high-S0s)/S0_high), color="blue", linewidth=lw)
    plt.loglog(n,10**N1(np.log10(n),S_popt1),color="orange",linestyle="--",linewidth=lw)
    #plt.loglog(n[:-1], abs((r0s[-1]-r0s[:-1])/r0s[-1]), linewidth=lw)
    plt.xlabel("Number of Cells")
    plt.ylabel("$\Sigma_{inner}$ Relative Difference")
    #plt.title("$\Sigma_{inner}$ Relative Difference vs. Number of Cells")
    ax = plt.gca()
    ax.set_box_aspect(1)
    plt.tight_layout()
    plt.show()
    
    return abs((S0s[-1]-S0s[:-1])/S0s[-1])



# Sigma Convergence Function

def ConvergenceMdotMax(Val_outer, Ms, Mp, Tp, Rp, r_in, r_out, alpha, n, Plot=True):
    
    Mdot_rel_diffs = []
    
    for i in range(len(n)):
        Mdot_max, Mdot_max_ana = MARC.Max_Acc_Rate(Val_outer, Ms, Mp, Tp, Rp, r_in, r_out, alpha, n[i],Const_Sigma=True)
        _, Sigma, _ = MARC.Max_Acc_Rate(Val_outer, Ms, Mp, Tp, Rp, r_in, r_out, alpha, n[i])
        Mdot_rel_diff = abs((Mdot_max-Mdot_max_ana)/Mdot_max_ana)
        Mdot_rel_diffs.append(Mdot_rel_diff)
    
    M_popt1, M_pcov1 = curve_fit(N1, np.log10(n), np.log10(Mdot_rel_diffs))

    Figsize, lw, Ticksize, Labelsize = [7, 2, 20, 20]
    
    plt.rc("figure", figsize=(Figsize,Figsize))
    plt.rc("xtick", labelsize=Ticksize)
    plt.rc("ytick", labelsize=Ticksize)
    plt.rc("axes", labelsize=Labelsize)
    
    # Convergence Plots
    plt.figure()
    plt.loglog(n, Mdot_rel_diffs, color="blue", linewidth=lw)
    plt.loglog(n,10**N1(np.log10(n),M_popt1),color="orange",linestyle="--",linewidth=lw)
    plt.xlabel("Number of Cells")
    plt.ylabel("$\dot{M}_{lim}$ Relative Difference")
    #plt.title("$\dot{M}_{max}$ Relative Difference vs. Number of Cells")
    ax = plt.gca()
    ax.set_box_aspect(1)
    plt.tight_layout()
    plt.show()
    
    return Mdot_rel_diffs



# Luminosity Finder IsoO Convergence Function

def ConvergenceLumIsoO(mc,me,rp,pd,kd,Td,L_guess,acc,n,m_grid_ini):
    
    M_E = 5.972e27                      # Earth mass [g]
    pc = 3.2                            # Core density [g/cm^3]
    rc = (3*mc*M_E/(4*np.pi*pc))**(1/3) # Core Radius [cm]
    
    Ls, rs, n_is = LFIO.L_Finder(mc,me,rp,pd,kd,Td,L_guess,acc,n,m_grid_ini, Plot=False)[2:]
    L_mean = np.mean(Ls[-4:])
    L_rel_diff = abs((np.array(Ls)-L_mean)/L_mean)
    r_rel_diff = abs((np.array(rs)-rc)/rc)
    
    Figsize, lw, Ticksize, Labelsize = [7, 2, 20, 20]
    
    plt.rc("figure", figsize=(Figsize,Figsize))
    plt.rc("xtick", labelsize=Ticksize)
    plt.rc("ytick", labelsize=Ticksize)
    plt.rc("axes", labelsize=Labelsize)
    
    plt.figure()
    plt.plot(n_is, L_rel_diff, linewidth=lw, color="blue")
    plt.yscale("log")
    plt.xlabel("Iteration Number")
    plt.ylabel("L Relative Difference")
    #plt.title("L Relative Difference Iteration")
    ax = plt.gca()
    ax.set_box_aspect(1)
    plt.tight_layout()
    plt.show()
    
    plt.figure()
    plt.plot(n_is, r_rel_diff, linewidth=lw, color="blue")
    plt.yscale("log")
    plt.xlabel("Iteration Number")
    plt.ylabel("$\epsilon_{r}$")
    #plt.title("$r_{core}$ Relative Difference Iteration")
    ax = plt.gca()
    ax.set_box_aspect(1)
    plt.tight_layout()
    plt.show()
    
    return L_rel_diff



# Luminosity Finder FreedmanO Convergence Function

def ConvergenceLumFreO(mc,me,rp,pd,kd,Td,L_guess,Z,acc_k,rc,acc,n,m_grid_ini):
    
    Ls, rs, n_is = LFFO.L_Finder(mc,me,rp,pd,kd,Td,L_guess,Z,acc_k,rc,acc,n,m_grid_ini, Plot=False)[3:]
    L_mean = np.mean(Ls[-3:])
    L_rel_diff = abs((Ls-L_mean)/L_mean)
    r_rel_diff = abs((np.array(rs)-rc)/rc)
    
    Figsize, lw, Ticksize, Labelsize = [7, 2, 20, 20]
    
    plt.rc("figure", figsize=(Figsize,Figsize))
    plt.rc("xtick", labelsize=Ticksize)
    plt.rc("ytick", labelsize=Ticksize)
    plt.rc("axes", labelsize=Labelsize)
    
    plt.figure()
    plt.plot(n_is, L_rel_diff, linewidth=lw)
    plt.yscale("log")
    plt.xlabel("Iteration Number")
    plt.ylabel("L Relative Difference")
    #plt.title("L Relative Difference Iteration")
    ax = plt.gca()
    ax.set_box_aspect(1)
    plt.tight_layout()
    plt.show()
    
    plt.figure()
    plt.plot(n_is, r_rel_diff, linewidth=lw)
    plt.yscale("log")
    plt.xlabel("Iteration Number")
    plt.ylabel("r Relative Difference")
    #plt.title("r Relative Difference Iteration")
    ax = plt.gca()
    ax.set_box_aspect(1)
    plt.tight_layout()
    plt.show()
    
    return L_rel_diff



# Luminosity Finder ComO Convergence Function

def ConvergenceLumComO(mc,me,rp,pd,kd,Td,L_guess,ad,d_to_g,acc_k,rc,acc,n,m_grid_ini):
    
    Ls, rs, n_is = LFCO.L_Finder(mc,me,rp,pd,kd,Td,L_guess,ad,d_to_g,acc_k,rc,acc,n,m_grid_ini, Plot=False)[3:]
    L_mean = np.mean(Ls[-3:])
    L_rel_diff = abs((Ls-L_mean)/L_mean)
    r_rel_diff = abs((np.array(rs)-rc)/rc)
    
    Figsize, lw, Ticksize, Labelsize = [7, 2, 20, 20]
    
    plt.rc("figure", figsize=(Figsize,Figsize))
    plt.rc("xtick", labelsize=Ticksize)
    plt.rc("ytick", labelsize=Ticksize)
    plt.rc("axes", labelsize=Labelsize)
    
    plt.figure()
    plt.plot(n_is, L_rel_diff, linewidth=lw)
    plt.yscale("log")
    plt.xlabel("Iteration Number")
    plt.ylabel("L Relative Difference")
    #plt.title("L Relative Difference Iteration")
    ax = plt.gca()
    ax.set_box_aspect(1)
    plt.tight_layout()
    plt.show()
    
    plt.figure()
    plt.plot(n_is, r_rel_diff, linewidth=lw)
    plt.yscale("log")
    plt.xlabel("Iteration Number")
    plt.ylabel("r Relative Difference")
    #plt.title("r Relative Difference Iteration")
    ax = plt.gca()
    ax.set_box_aspect(1)
    plt.tight_layout()
    plt.show()
    
    return L_rel_diff
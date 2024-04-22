#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 16:46:31 2022

@author: ollie
"""

import numpy as np
import IsoO_Model as IOM
import matplotlib.pyplot as plt





def L_Finder(mc,me,rp,pd,kd,Td,L_guess,acc,n,m_grid_ini, Plot=True):
    
    M_E, R_J = 5.972e27, 6.9911e9           #Earth mass [g], Jupiter radius [cm]
    pc = 3.2                                # Core density [g/cm^3]
    rc = (3*mc*M_E/(4*np.pi*pc))**(1/3)
    
    L = L_guess
    L_step = 1/2
    e_r = 1.
    n_i = 0
    
    Ls, rs, n_is, Plus_Minus = [], [], [], []
    
    while abs(e_r) > acc:
        
        n_i += 1
        n_is.append(n_i)
        
        IsoO = IOM.Iso_Opacity_Giant_Planet(mc,me,rp,pd,kd,Td,L,n,grid_given=True,m_grid=m_grid_ini)
        var = IsoO.Planet_Formation("Both")
        r_inner = var[1][1]
        
        Ls.append(L)
        rs.append(r_inner)
        
        e_r = r_inner/rc - 1
        
        if abs(e_r) >  acc:
            if e_r > 0:
                Plus_Minus.append("+")
                if n_i > 1:
                    if Plus_Minus[-2] == "+":
                        L_step *= 2
                L *= (1 + L_step)
            elif e_r < 0:
                Plus_Minus.append("-")
                if n_i > 1:
                    if Plus_Minus[-2] == "-":
                        L_step *= 2
                L *= (1 - L_step)
        
        elif abs(e_r) <= acc:
            L = L
        
        L_step /= 2
    
    Ls = np.array(Ls)
    rs = np.array(rs)
    n_is = np.array(n_is)
    
    if Plot == True:
        Figsize, lw, Ticksize, Labelsize, Legendfontsize = [7, 2, 20, 20, 15]
        
        plt.rc("figure", figsize=(Figsize,Figsize))
        plt.rc("xtick", labelsize=Ticksize)
        plt.rc("ytick", labelsize=Ticksize)
        plt.rc("axes", labelsize=Labelsize)
        
        plt.figure()
        plt.plot(n_is, Ls/1e-6, linewidth=lw, color="blue")
        plt.xlabel("Iteration Number")
        plt.ylabel("Luminosity / 10$^{-6}$ L$_{\odot}$")
        #plt.title("L Iteration")
        ax = plt.gca()
        ax.set_box_aspect(1)
        plt.tight_layout()
        plt.show()
        
        plt.figure()
        plt.plot(n_is, rs/R_J, linewidth=lw, color="blue")
        plt.axhline(rc/R_J, linestyle="--", color="black")
        plt.xlabel("Iteration Number")
        plt.ylabel("$r_{inner}$ / $R_{J}$")
        #plt.title("r Iteration")
        ax = plt.gca()
        ax.set_box_aspect(1)
        ax = plt.gca()
        ax.set_yticks([0,0.1,0.2,rc/R_J,0.3,0.4])
        ax.set_yticklabels([0,0.1,0.2,"$r_{core}$",0.3,0.4])
        plt.tight_layout()
        plt.show()
        
        fig, ax1 = plt.subplots()
        color1 = 'red'
        color2 = 'orange'
        
        ax1.set_xlabel("Iteration Number")
        ax1.set_ylabel("Luminosity / 10$^{-6}$ $L_{\odot}$", color=color1)
        ax1.plot(n_is, Ls/1e-6, color=color1, linewidth=lw)
        ax1.tick_params(axis='y')
        
        ax2 = ax1.twinx()
        ax2.set_ylabel("$r_{inner}$ / $R_{J}$", color=color2)
        ax2.plot(n_is, rs/R_J, color=color2, linewidth=lw)
        ax2.axhline(rc/R_J, linestyle="--", color="black", linewidth=lw)
        ax2.tick_params(axis='y')
        fig.tight_layout()
    
    return L, var, Ls, rs, n_is
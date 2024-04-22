#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 16:46:31 2022

@author: ollie
"""

import numpy as np
import GPF_Freedman_Opacity as GPFFO
import matplotlib.pyplot as plt





def L_Finder(mc,me,rp,pd,kd,Td,L_guess,acc,Z,acc_k,n,m_grid_ini, Plot=True):
    
    M_E = 5.972e27                  # Earth mass [g]
    pc = 3.2                        # Core density [g/cm^3]
    rc = (3*mc*M_E/(4*np.pi*pc))**(1/3)
    
    L = L_guess
    L_step = 1/2
    e_r = 1.
    n_i = 0
    
    Ls, rs, n_is, Plus_Minus, breaks = [], [], [], [], []
    
    while abs(e_r) > acc:
        
        n_i += 1
        n_is.append(n_i)
        
        FreedmanO = GPFFO.Freedman_Opacity_Giant_Planet(mc,me,rp,pd,kd,Td,L,Z,acc_k,n,grid_given=True,m_grid=m_grid_ini)
        var = FreedmanO.Planet_Formation()
        r_inner = var[1][1]
        
        Ls.append(L)
        rs.append(r_inner)
        
        breaks.append(FreedmanO._broken)
        
        if breaks[-1] == True:
            Plus_Minus.append("-")
            L *= (1 - L_step)
            continue
        
        e_r = r_inner/rc - 1
        
        if abs(e_r) >  acc:
            if e_r > 0:
                Plus_Minus.append("+")
                if n_i > 1:
                    if Plus_Minus[-2] == "+":
                        L_step *= 4
                L *= (1 + L_step)
            elif e_r < 0:
                Plus_Minus.append("-")
                if n_i > 1:
                    if Plus_Minus[-2] == "-":
                        L_step *= 4
                L *= (1 - L_step)
        
        elif abs(e_r) <= acc:
            L = L
        
        #print(L, r_inner/rc -1)
        
        L_step /= 4
    
    
    if Plot == True:
        plt.figure()
        plt.plot(n_is, Ls, "x")
        plt.xlabel("Iteration Number")
        plt.ylabel("L")
        plt.title("L Iteration")
        plt.show()
        
        plt.figure()
        plt.plot(n_is, rs, "x")
        plt.xlabel("Iteration Number")
        plt.ylabel("r")
        plt.title("r Iteration")
        plt.show()
    
    return L, var, FreedmanO, Ls, rs, n_is
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 16:46:31 2022

@author: ollie
"""

import GPF_Combined_Opacity as GPFCO
import matplotlib.pyplot as plt





def L_Finder(mc,me,rp,pd,kd,Td,L_guess,ad,d_to_g,acc_k,rc,acc,n,m_grid_ini, Plot=True):
    
    L = L_guess
    L_step = 1/2
    e_r = 1.
    n_i = 0
    
    Ls, rs, n_is, Plus_Minus, breaks = [], [], [], [], []
    
    while abs(e_r) > acc:
        
        n_i += 1
        n_is.append(n_i)
        
        CombinedO = GPFCO.Freedman_Opacity_Giant_Planet(mc,me,rp,pd,kd,Td,L,ad,d_to_g,acc_k,n,grid_given=True,m_grid=m_grid_ini)
        var = CombinedO.Planet_Formation()
        r_inner = var[1][1]
        
        Ls.append(L)
        rs.append(r_inner)
        
        breaks.append(CombinedO._broken)
        
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
    
    return L, var, CombinedO, Ls, rs, n_is
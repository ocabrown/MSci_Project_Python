#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 19 16:52:33 2022

@author: ollie
"""




import IsoT_Model as ITM
import IsoD_Model as IDM
import IsoO_Model as IOM
import GPF_Movshovitz_Opacity_Forced as GPFMOF
import GPF_Movshovitz_Opacity as GPFMO
import GPF_Freedman_Opacity as GPFFO
import GPF_Combined_Opacity as GPFCO



def Mass_Grid_Solver_Iso_Temperature(mc,me,rp,pd,Td,n,acc_m,p_ana_c_ini):
    
    e_m0 = 1.
    e_m1 = 1.
    p_ana_c_guess = p_ana_c_ini
    
    while (e_m0 > acc_m) or (e_m1 > acc_m):
        
        IsoT = ITM.Iso_Temperature_Giant_Planet(mc,me,rp,pd,Td,n,p_ana_c=p_ana_c_guess)
        IsoT.Planet_Formation()
        p_ana_c_guess_new = IsoT.DensityFit()
        
        e_m0 = abs(p_ana_c_guess_new[0]/p_ana_c_guess[0] - 1)
        e_m1 = abs(p_ana_c_guess_new[1]/p_ana_c_guess[1] - 1)
        
        p_ana_c_guess = p_ana_c_guess_new
    
    m_grid = IsoT._var[0]
    
    return m_grid, p_ana_c_guess_new


"""
def Mass_Grid_Solver_Iso_Density(mc,me,rp,pd,kd,Td,lp,n,acc_m,p_ana_c_ini):
    
    e_m0 = 1.
    e_m1 = 1.
    p_ana_c_guess = p_ana_c_ini
    
    while (e_m0 > acc_m) or (e_m1 > acc_m):
        
        IsoD = IDM.Iso_Density_Giant_Planet(mc,me,rp,pd,kd,Td,lp,n,p_ana_c=p_ana_c_guess)
        IsoD.Planet_Formation()
        p_ana_c_guess_new = IsoD.DensityFit()
        
        e_m0 = abs(p_ana_c_guess_new[0]/p_ana_c_guess[0] - 1)
        e_m1 = abs(p_ana_c_guess_new[1]/p_ana_c_guess[1] - 1)
        
        p_ana_c_guess = p_ana_c_guess_new
    
    m_grid = IsoD._var[0]
    
    return m_grid, p_ana_c_guess_new
"""


def Mass_Grid_Solver_Iso_Opacity(mc,me,rp,pd,kd,Td,lp,n,acc_m,p_ana_c_ini):
    
    e_m0 = 1.
    e_m1 = 1.
    p_ana_c_guess = p_ana_c_ini
    
    while (e_m0 > acc_m) or (e_m1 > acc_m):
        
        IsoO = IOM.Iso_Opacity_Giant_Planet(mc,me,rp,pd,kd,Td,lp,n,p_ana_c=p_ana_c_guess)
        IsoO.Planet_Formation()
        p_ana_c_guess_new = IsoO.DensityFit()
        
        e_m0 = abs(p_ana_c_guess_new[0]/p_ana_c_guess[0] - 1)
        e_m1 = abs(p_ana_c_guess_new[1]/p_ana_c_guess[1] - 1)
        
        p_ana_c_guess = p_ana_c_guess_new
    
    m_grid = IsoO._var[0]
    
    return m_grid, p_ana_c_guess_new



def Mass_Grid_Solver_Movshovitz_Opacity_Forced(mc,me,rp,pd,kd,Td,lp,n,acc_m,p_ana_c_ini,Data):
    
    e_m0 = 1.
    e_m1 = 1.
    p_ana_c_guess = p_ana_c_ini
    
    while (e_m0 > acc_m) or (e_m1 > acc_m):
        
        MovOF = GPFMOF.Movshovitz_Opacity_Forced_Giant_Planet(mc,me,rp,pd,kd,Td,lp,n,p_ana_c=p_ana_c_guess)
        MovOF.Planet_Formation(Data)
        p_ana_c_guess_new = MovOF.DensityFit()
        
        e_m0 = abs(p_ana_c_guess_new[0]/p_ana_c_guess[0] - 1)
        e_m1 = abs(p_ana_c_guess_new[1]/p_ana_c_guess[1] - 1)
        
        p_ana_c_guess = p_ana_c_guess_new
    
    m_grid = MovOF._var[0]
    
    return m_grid, p_ana_c_guess_new



def Mass_Grid_Solver_Movshovitz_Opacity(mc,me,rp,pd,kd,Td,lp,n,acc_m,p_ana_c_ini,Data):
    
    e_m0 = 1.
    e_m1 = 1.
    p_ana_c_guess = p_ana_c_ini
    
    while (e_m0 > acc_m) or (e_m1 > acc_m):
        
        MovO = GPFMO.Movshovitz_Opacity_Giant_Planet(mc,me,rp,pd,kd,Td,lp,n,p_ana_c=p_ana_c_guess)
        MovO.Planet_Formation(Data)
        p_ana_c_guess_new = MovO.DensityFit()
        
        e_m0 = abs(p_ana_c_guess_new[0]/p_ana_c_guess[0] - 1)
        e_m1 = abs(p_ana_c_guess_new[1]/p_ana_c_guess[1] - 1)
        
        p_ana_c_guess = p_ana_c_guess_new
    
    m_grid = MovO._var[0]
    
    return m_grid, p_ana_c_guess_new



def Mass_Grid_Solver_Freedman_Opacity(mc,me,rp,pd,kd,Td,lp,met,acc_k,n,acc_m,p_ana_c_ini):
    
    e_m0 = 1.
    e_m1 = 1.
    p_ana_c_guess = p_ana_c_ini
    
    while (e_m0 > acc_m) or (e_m1 > acc_m):
        
        FreO = GPFFO.Freedman_Opacity_Giant_Planet(mc,me,rp,pd,kd,Td,lp,met,acc_k,n,p_ana_c=p_ana_c_guess)
        FreO.Planet_Formation()
        p_ana_c_guess_new = FreO.DensityFit()
        
        e_m0 = abs(p_ana_c_guess_new[0]/p_ana_c_guess[0] - 1)
        e_m1 = abs(p_ana_c_guess_new[1]/p_ana_c_guess[1] - 1)
        
        p_ana_c_guess = p_ana_c_guess_new
    
    m_grid = FreO._var[0]
    
    return m_grid, p_ana_c_guess_new



def Mass_Grid_Solver_Combined_Opacity(mc,me,rp,pd,kd,Td,lp,ad,d_to_g,acc_k,n,acc_m,p_ana_c_ini):
    
    e_m0 = 1.
    e_m1 = 1.
    p_ana_c_guess = p_ana_c_ini
    
    while (e_m0 > acc_m) or (e_m1 > acc_m):
        
        ComO = GPFCO.Combined_Opacity_Giant_Planet(mc,me,rp,pd,kd,Td,lp,ad,d_to_g,acc_k,n,p_ana_c=p_ana_c_guess)
        ComO.Planet_Formation()
        p_ana_c_guess_new = ComO.DensityFit()
        
        e_m0 = abs(p_ana_c_guess_new[0]/p_ana_c_guess[0] - 1)
        e_m1 = abs(p_ana_c_guess_new[1]/p_ana_c_guess[1] - 1)
        
        p_ana_c_guess = p_ana_c_guess_new
        print(p_ana_c_guess_new)
    
    m_grid = ComO._var[0]
    
    return m_grid, p_ana_c_guess_new
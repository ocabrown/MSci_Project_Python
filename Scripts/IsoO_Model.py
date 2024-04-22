#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 19 16:52:33 2022

@author: ollie
"""

import numpy as np
import scipy.constants as spc
import scipy.optimize as so





class Iso_Opacity_Giant_Planet:
    
    
    
    def __init__(self,mc,me,rp,pd,kd,Td,lp,n,grid_given=False,m_grid=None,p_ana_c=None):
        
        # Constants:
        num_p = int(n)                  # Number of points along m
        num_g = 1                       # Number of ghost cells on each side
        
        AU = spc.au * 1e2               # Astronomical unit [cm]
        G = spc.G * 1e3                 # Gravitational constant
        kB = spc.k * 1e7                # Boltzmann constant
        mu = 2.3                        # Relative mass
        R = 8.31e7                      # Gas constant [erg/mol/K]
        R /= mu
        mp = spc.m_p * 1e3              # Proton mass [g]
        R_J = 6.9911e9                  # Jupiter radius [cm]
        L_Sun = 3.846e33                # Solar luminosity [erg/s]
        M_E = 5.972e27                  # Earth mass [g]
        M_S = 1.989e33                  # Solar mass [g]
        
        cs2 = R * Td            # Sound speed [cm/s]
        
        self._lo = num_g
        self._hi = num_g + num_p - 1
        
        core_mass_max = mc * M_E        # Planet's core mass [M_E]
        env_mass_max = me * M_E         # Planet's envelope mass [M_E]
        Mp = core_mass_max + env_mass_max
        R_H = (5.2*AU) * (Mp/(3*M_S))**(1/3)
        #rad_hi = (G*Mp/((cs2/1.)+((G*Mp)/(0.25*R_H))))
        #print(rad_hi/R_J)
        rad_hi = rp * R_J               # Planet radius [cm]
        den_hi = pd                     # Disc density [g/cm^3]
        opa_hi = kd                     # Disc opacity [cm^2/g]
        tem_hi = Td                     # Disc temperature [K]
        pre_hi = pd * R * Td            # Disc pressure [g/cm/s^2]
        lum = lp * L_Sun                # Luminosity [L_Sun]
        
        
        # Making empty arrays for other variables with boundary conditions
        empty = np.zeros((num_p + (2*num_g)))
        r, p, k, T, P, s = np.array([empty,empty,empty,empty,empty,empty])
        r[self._hi] = rad_hi
        p[self._hi] = den_hi
        k[self._hi] = opa_hi
        T[self._hi] = tem_hi
        P[self._hi] = pre_hi
        
        
        if grid_given == False:
            # Calculating envelope mass array so r has Uniform logarithmic spacing
            env_mass_min = 1e20         # [g]
            
            def m_r_func(r):
                a = p_ana_c[0]
                b = p_ana_c[1]
                integral = -(4*np.pi*((a*r)**2 + 2*a*r + 2) * np.e**(-(a*r+b)))/a**3
                C = env_mass_max - (-(4*np.pi*((a*rad_hi)**2 + 2*a*rad_hi + 2) * np.e**(-(a*rad_hi+b)))/a**3)
                return integral + C
            
            def m_min(r):
                return env_mass_min - m_r_func(r)
            
            
            rad_lo_guess = np.linspace(1e8, 1e12, 10000)
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
        
        else:
            m = m_grid
        
        """
        # Uniform logarithmic spacing for me (envelope mass)
        env_mass_min = 1e20
        m1 = np.geomspace(env_mass_min, 0.9*env_mass_max, num=int(0.35*num_p))
        m2 = np.geomspace((0.9*env_mass_max), env_mass_max, num=int((0.65*num_p)+1))
        m2 = m2[1:]
        m = np.concatenate((m1,m2))
        for i in range(num_g):
            m = np.insert(m, 0,      0.)
            m = np.insert(m, len(m), 0.)
        """
        
        self._const = [core_mass_max, lum, R]
        self._var_ini = [m, r, p, k, T, P, s]
    
    
    
    def Planet_Formation(self, Rad_or_Conv="Both"):
        
        # Constants:
        G = spc.G * 1e3                 # Gravitational constant
        a = 7.57e-15                    # erg/cm^3/K^4
        c = spc.c * 1e2                 # Light speed [cm/s]
        gamma = 7/5                     # Adiabatic constant
        M_core = self._const[0]
        l = self._const[1]
        R = self._const[2]              # Gas constant [erg/mol/K]
        
        Grad_Adia = (gamma-1)/gamma     # Adiabatic temperature gradient term
        
        var_new = self._var_ini.copy()
        m, r, p, k, T, P, s = var_new
        dP_dr_num = s.copy()
        dP_dr_ana = s.copy()
        
        self._conv_or_rad = []
        self._boundary_i = []
        
        for i in range(self._lo+1, self._hi+1):
            
            r[-(i+1)] = (r[-i]**3 + 3*(m[-(i+1)]-m[-i])/(4*np.pi*p[-i]))**(1/3)
            P[-(i+1)] = P[-i] - G*(m[-(i+1)]+M_core)*(m[-(i+1)]-m[-i])/(4*np.pi*r[-(i+1)]**4)
            k[-(i+1)] = k[-i]
            
            if Rad_or_Conv == "Purely Radiative":
                T[-(i+1)] = (T[-i]**4 - 3*k[-(i+1)]*l*(m[-(i+1)]-m[-i])/(16*a*c*(np.pi**2)*r[-(i+1)]**4))**(1/4)
            
            elif Rad_or_Conv == "Purely Convective":
                T[-(i+1)] = T[-i] * ((P[-(i+1)]/P[-i])**((gamma-1)/gamma))
            
            elif Rad_or_Conv == "Both":
                Grad_Rad = 3*k[-(i+1)]*l*P[-i]/(16*np.pi*a*c*G*(m[-(i+1)]+M_core)*(T[-i]**4))
                if abs(Grad_Rad) >= abs(Grad_Adia):         # Convective
                    self._conv_or_rad.append(0)
                    T[-(i+1)] = T[-i] * ((P[-(i+1)]/P[-i])**((gamma-1)/gamma))
                
                else:                                       # Radiative
                    self._conv_or_rad.append(1)
                    T[-(i+1)] = (T[-i]**4 - 3*k[-(i+1)]*l*(m[-(i+1)]-m[-i])/(16*a*c*(np.pi**2)*r[-(i+1)]**4))**(1/4)
                
                if i > 2:
                    if self._conv_or_rad[-1] - self._conv_or_rad[-2] == 1:
                        self._boundary_i.append([i,1])
                    elif self._conv_or_rad[-1] - self._conv_or_rad[-2] == -1:
                        self._boundary_i.append([i,0])
            
            p[-(i+1)] = P[-(i+1)]/(R*T[-(i+1)])
            s[-(i+1)] = T[-(i+1)]/(P[-(i+1)]**((gamma-1)/gamma))
            dP_dr_num[-(i+1)] = (P[-(i+1)]-P[-i])/(r[-(i+1)]-r[-i])
            dP_dr_ana[-(i+1)] = -G*(m[-(i+1)]+M_core)*p[-(i+1)]/(r[-(i+1)]**2)
        
        self._varL2 = [m, r, p, k, T, P, s, dP_dr_num, dP_dr_ana]
        self._var = [m, r, p, k, T, P, s]
        
        return self._var
    
    
    
    def L2norm(self):
        
        var_conv = self._varL2.copy()
        m, r, p, k, T, P, s, dP_dr_num, dP_dr_ana = [var[self._lo:self._hi] for var in var_conv]
        
        min_pos = np.argmin(dP_dr_num)
        min_dP_dr = dP_dr_num[min_pos]
        
        L2norm_dP_dr = np.sqrt(sum(((dP_dr_ana - dP_dr_num)/dP_dr_ana)**2)/len(dP_dr_num))
        L2norms_dP_dr = ((dP_dr_ana - dP_dr_num)/dP_dr_ana)**2
        
        return L2norm_dP_dr, min_dP_dr, L2norms_dP_dr
    
    
    
    def DensityFit(self):
        
        var_plot = self._var.copy()
        m, r, p, k, T, P, s = [var[self._lo:self._hi] for var in var_plot]
        
        def func(x,a,b):
            return -(a*x + b)
        popt, pcov = so.curve_fit(func,r,np.log(p))
        
        return popt
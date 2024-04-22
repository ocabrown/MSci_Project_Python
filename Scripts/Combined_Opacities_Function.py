#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 11 18:44:23 2023

@author: ollie
"""

import numpy as np
import scipy.constants as spc


###############################################################################
# Opacity-fit from Freedman et al. 2014
#      Input:
#           P in cgs   e.g. 1e6 is 1 bar
#           T in cgs   e.g. 100 is 100 K
#           Z in log([M/H] / [M/H]_solar), e.g. 0 is solar
#      Output:
#           Rosseland opacity of the dust/gas mixture in cgs, i.e. cm^2/g
###############################################################################

def Freedman_Opacity(P, T, Z): 

    # Values from Table 2 of paper:
    c  = [0., 10.602, 2.882, 6.09e-15, 2.954, -2.526, 0.843, -5.490]
    if T < 800:
        d = [-14.051, 3.055, 0.024, 1.877, -0.4451766673621382, 0.8321]
    else:
        d = [82.241, -55.456, 8.753389704734504, 0.7048, -0.0414, 0.8321]
    
    logT = np.log10(T)       # Log of temperature
    logP = np.log10(P)       # Log of pressure
    
    # Equations (4), (5) and (3) from paper:
    logkloP = c[1] * np.arctan(logT - c[2]) - (c[3] / (logP + c[4])) * np.exp((logT - c[5])*(logT - c[5]))+ c[6] * Z + c[7] 
    logkhiP = d[0] + d[1]*logT + d[2]*logT*logT + (d[3] + d[4]*logT)*(logP) + d[5] * Z * (0.5 + np.arctan(5*logT - 12.5)/np.pi) 
    kgas = 10.**logkloP + 10.**logkhiP
    
    return kgas



###############################################################################
#
#
# Mean opacity - Rosseland
#
#
###############################################################################

h = spc.h * 1e7
c = spc.c * 1e2
k = spc.k * 1e7
planck_a = 2. * h * (c**2)
planck_b = h * c / k

###############################################################################
# Planck B(T) = (2hc^2/l^5) / (e^(hc/lkT) - 1)
# dB(T)/dT = (2h^2c^3/l^6kT^2) * (e^(hc/lkT)/((e^(hc/lkT)-1)^2))
# wav: wavelength in cm
# T: temperature in K
# returns: Energy dens in erg cm^-3 s^-1 sr^-1 K^-1

def Planck_Deriv(l, T):
    dB_dT = (2.*(h**2)*(c**3)/((l**6)*k*(T**2))) * (1. / ((np.exp(h*c/(l*k*T))-1.) * (1.-np.exp(-h*c/(l*k*T)))))
    return dB_dT



###############################################################################
# Dust effective absorption function
# x: wavelength in cm
# returns: Exctinction efficiency [1], see e.g. Mordasini 2014b (Analytical dust model paper) and refs therein

def Q(x):
    if x < 1e-3:
        return 2.
    elif x < 0.5:
        return 2.+4.*x
    elif x > 8./3.:
        return 0.3 * x**(-1.)
    else:
        return 0.8 * x**(-2.)



###############################################################################
# Wrapper function
# l: wavelength in cm
# a: dust size  in cm
# output: opacity in cm^2 per gram of dust at all given x

def dust_opa(ls, a):
    xs = ls/(2.*np.pi*a)
    return np.array([Q(x)/a for x in xs ])
"""
####################################################################################################
# Dust effective absorption function
# x: 2 pi a / l
# returns: Exctinction efficiency [1], see e.g. Mordasini 2014b eqn 54 (Analytical dust model paper) and refs therein

def Q(x):
    
    if x < 0.375:
        return 0.3 * x
    elif 0.375 <= x < 2.188:
        return 0.8 * (x**2)
    elif 2.188 <= x < 1000:
        return 2. + (4./x)
    else:
        return 2.



####################################################################################################
# Wrapper function
# l: wavelength in cm
# a: dust size  in cm
# output: opacity in cm^2/g of dust at all given x

def dust_opa(l, a):
    xs = 2. * np.pi * a / l
    return np.array([Q(x)/a for x in xs])
"""


###############################################################################
# On the fly Rosseland mean opacity for dust of a given size a and temperature field T
# T: temperature [K]
# a: particle size [cm]
# returns: Rosseland mean opacity, cm^2/g

def Dust_Ross(T,a):
    ls  = 10.**np.linspace(-4,4,1000) * 1e-4        # mum to cm
    dls= np.diff(np.hstack((0, ls)))

    k_dust = dust_opa(ls,a)
    norm = np.sum([Planck_Deriv(l,T)*dls[li] for li,l in enumerate(ls)])
    dust_norm = np.sum([(1/k_dust[li])*Planck_Deriv(l,T)*dls[li] for li,l in enumerate(ls)])
    
    if T < 1500.:
        return norm/dust_norm
    else:
        return 0.                                   # Due to evaporation





###############################################################################
#
#
# Combined Rosseland mean opacity for dust and gas under a Rosseland addition rule assumption
#
#
###############################################################################

###############################################################################
# P: gas pressure [dyne/cm^2], ref value: 1e6=1bar
# T: temperature [K]
# a: particle size [cm]
# dust_to_gas: mass ratio of dust to gas, sets metallicity for gas opacities as well
# returns: cm^2/g

def Combined_Opacity(P,T,a,dust_to_gas):
    
    gas_opacity = Freedman_Opacity(P, T, np.log10(dust_to_gas/1e-2))
    dust_opacity = Dust_Ross(T,a) * dust_to_gas
    opacity = gas_opacity + dust_opacity
    
    return opacity





"""
T = 10.**np.linspace(1,4.1,100)
P_val = 1e6
P = np.full((500), P_val)
k = []
for i in range(len(T)):
    k.append(co.combined_opacity(P[i],T[i],1.5e-4,1e-2))
k = np.array(k)
plt.figure()
plt.loglog(T,k)
plt.xlabel("Temperature / K")
plt.ylabel("Opacity")
plt.show()

print(k)
"""


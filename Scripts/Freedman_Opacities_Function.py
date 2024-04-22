#
# Opacity-fit from Freedman et al. 2014
#      Input:
#           P in cgs   e.g. 1e6 is 1 bar
#           T in cgs   e.g. 100 is 100 K
#           Z in log([M/H] / [M/H]_solar), e.g. 0 is solar
#      Output:
#           Rosseland opacity of the dust/gas mixture in cgs, i.e. cm^2/g
#

import numpy as np



def Freedman_Opacity(Z, P, T): 
    
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
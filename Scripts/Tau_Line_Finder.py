#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 15:03:53 2023

@author: ollie
"""

import numpy as np



def Tau_Line_Finder(r, p, k, tau):
    
    height = p*k
    width = r
    area_list = []
    
    for i in range(len(width)+1): #integrating for optical depth, stopping when tau=1
        area_val = (width[-(i+1)]-width[-(i+2)])*height[-(i+1)]
        area_list.append(area_val)
        area = np.sum(area_list)
        if area > tau:
            tau_i = i
            break
    
    return tau_i
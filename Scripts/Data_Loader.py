#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 19:56:50 2022

@author: ollie
"""

import numpy as np



# Comparison Data Extractor

def Data_Extractor(Data, Type):
    
    X_Data = 10**Data[:,0]
    Y_Data = Data[:,1]  
    
    if Type == "Density":
        Y_Data = 10**Y_Data[:]
    elif Type == "Opacity":
        Y_Data = 10**Y_Data[:]
    elif Type == "Temperature":
        Y_Data = Y_Data
    
    return X_Data, Y_Data



def Data_Preparer(Filenames):
    
    # Import Data
    DensityData = np.loadtxt("/Data/"+str(Filenames[0]), delimiter=",")
    OpacityData = np.loadtxt("/Data/"+str(Filenames[1]), delimiter=",")
    TemperatureData = np.loadtxt("/Data/"+str(Filenames[2]), delimiter=",")
    RD, D = Data_Extractor(DensityData, "Density")
    RO, O = Data_Extractor(OpacityData, "Opacity")
    RT, T = Data_Extractor(TemperatureData, "Temperature")
    
    return [[RD,D], [RO,O], [RT,T]]



def Data_Extractor_Sigma(Filename):
    
    SigmaData = np.loadtxt("/Data/"+str(Filename), delimiter=",")
    X_Data = SigmaData[:,0]
    Y_Data = SigmaData[:,1]
    
    return X_Data, Y_Data
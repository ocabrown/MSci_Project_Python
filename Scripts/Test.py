#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 19 22:00:46 2022

@author: ollie
"""

import numpy as np
import time

import Data_Loader as DL
import Mass_Grid_Solver as MGS
import Convergence_Tester as CT
import IsoT_Model as ITM
import IsoD_Model as IDM
import IsoO_Model as IOM
import GPF_Movshovitz_Opacity_Forced as GPFMOF
import GPF_Movshovitz_Opacity as GPFMO
import GPF_Freedman_Opacity as GPFFO
import GPF_Combined_Opacity as GPFCO
import Maximum_Accretion_Rate_Calculator as MARC
import L_Finder_IsoO as LFIO
import Evolver_IsoO as EIO
import L_Finder_FreO as LFFO
import Evolver_FreO as EFO
import L_Finder_ComO as LFCO
import Evolver_ComO as ECO
import Plotter





#%% Loading Comparison Data

Movshovitz_Figure_3 = DL.Data_Preparer(["Movshovitz Data/Figure 3 - Density.csv", "Movshovitz Data/Figure 3 - Opacity.csv", "Movshovitz Data/Figure 3 - Temperature.csv"])
Movshovitz_Figure_3_Opacity_Log = np.log(Movshovitz_Figure_3[1])
Movshovitz_Figure_4 = DL.Data_Preparer(["Movshovitz Data/Figure 4 - Density.csv", "Movshovitz Data/Figure 4 - Opacity.csv", "Movshovitz Data/Figure 4 - Temperature.csv"])
Movshovitz_Figure_4_Opacity_Log = np.log(Movshovitz_Figure_4[1])

Crida_Sigma_Data = DL.Data_Extractor_Sigma("Crida Data/Sigma (0.05, -5.5).csv")





#%%


#                       Constant Temperature Testing



#%% Values are just test values to show the code works (not necessarily realisitic)
# mc, me, rp, pd, Td, n, acc_m, p_ana_c_ini
val_IsoT = [10., 0.01, 41.4089, 1.47e-11, 118, 1e3, 1e-6, [4.11e-11, 1.49e1]] #[4.58e-11, 1.46e1], 36.3776
mc_IsoT, me_IsoT, rp_IsoT, pd_IsoT, Td_IsoT, n_IsoT, acc_m, p_ana_c_ini_IsoT = val_IsoT



#%% Calculating the Mass Grid

m_grid_IsoT, _ = MGS.Mass_Grid_Solver_Iso_Temperature(mc_IsoT,me_IsoT,rp_IsoT,pd_IsoT,Td_IsoT,n_IsoT,acc_m,p_ana_c_ini_IsoT)



#%% Running Simulation

IsoT1 = ITM.Iso_Temperature_Giant_Planet(mc_IsoT,me_IsoT,rp_IsoT,pd_IsoT,Td_IsoT,n_IsoT,grid_given=True,m_grid=m_grid_IsoT)
IsoT1_Sim_Data = IsoT1.Planet_Formation()
print("Max p/p_Disc =", max(IsoT1_Sim_Data[2])/pd_IsoT)

IsoT1.Plot("p_comp", "r")
IsoT1.Plot("dP/dr", "r")



#%% Extra Plots

IsoT1.Plot("m", "r")
IsoT1.Plot("P", "r")



#%% Viva plot

IsoT1.Plot("sharedx", "r")



#%% Convergence Test - can vary n_array to see the numerical precision

n_array_IsoT = np.linspace(1e3,1e4,10)
#n_array_IsoT = np.concatenate((np.linspace(50,1e2,10), np.linspace(1e2,1e3,10)))#, np.linspace(1e3,1e4,10)))
L2norms_p_IsoT, L2norms_dP_dr_upper_IsoT = CT.ConvergenceIsoT(mc_IsoT,me_IsoT,rp_IsoT,pd_IsoT,Td_IsoT,n_array_IsoT,acc_m,p_ana_c_ini_IsoT)



#%% L2norms full range 1e2 - 1e5

n_array_IsoT15 = np.concatenate((np.linspace(1e2,1e3,30), np.linspace(1e3,1e4,50), np.linspace(1e4,1e5,30)))
L2norms_p_IsoT15, L2norms_dP_dr_upper_IsoT15 = CT.ConvergenceIsoT(mc_IsoT,me_IsoT,rp_IsoT,pd_IsoT,Td_IsoT,n_array_IsoT15,acc_m,p_ana_c_ini_IsoT)



#%% Shows the increase is in fact numerical noise as L2norm per cell average around 1e-5 so for 1e7 cells L2norm is 1e2

n_IsoT2 = 1e6
IsoT2 = ITM.Iso_Temperature_Giant_Planet(mc_IsoT,me_IsoT,rp_IsoT,pd_IsoT,Td_IsoT,n_IsoT2,p_ana_c=p_ana_c_ini_IsoT)
IsoT2_Sim_Data = IsoT2.Planet_Formation()

L2norms_p2, L2norms_dP_dr2 = IsoT2.L2norm()[3:]
avr_L2norm_p2, max_L2norm_p2 = np.mean(L2norms_p2), max(L2norms_p2)
avr_L2norm_dP_dr2, max_L2norm_dP_dr2 = np.mean(L2norms_dP_dr2), max(L2norms_dP_dr2)
print("Density avr. and max. L2norm =", avr_L2norm_p2, "and", max_L2norm_p2)
print("Press. grad. avr. and max. L2norm =", avr_L2norm_dP_dr2, "and", max_L2norm_dP_dr2)





#%%


#                       Constant Density Testing



#%% Values are just test values to show the code works (not necessarily realisitic)
# mc, me, rp, pd, kd, Td, lp, n
#val_IsoD = [10., 0.0001, 30, 1e-10, 100, 118, 10, 1e3]
val_IsoD = [10., 0.0008, 36.3776, 7e-11, 1.5e-5, 118, 10000, 1e3] #38.9459
mc_IsoD, me_IsoD, rp_IsoD, pd_IsoD, kd_IsoD, Td_IsoD, lp_IsoD, n_IsoD = val_IsoD



#%% Running Simulation

IsoD1 = IDM.Iso_Density_Giant_Planet(mc_IsoD,me_IsoD,rp_IsoD,pd_IsoD,kd_IsoD,Td_IsoD,lp_IsoD,n_IsoD)
IsoD1_Sim_Data = IsoD1.Planet_Formation()
print("Max T/T_Disc =", max(IsoD1_Sim_Data[3])/Td_IsoD)

IsoD1.Plot("T_comp", "r")
IsoD1.Plot("dP/dr", "r")
#IsoD1.Plot("dP/dr_lower", "r")
#IsoD1.Plot("dP/dr_upper", "r")



#%% Extra Plots

IsoD1.Plot("m", "r")
IsoD1.Plot("P", "r")



#%% Viva Plot

IsoD1.Plot("sharedx", "r")



#%% Convergence Test - can vary n_array to see the numerical precision

n_array_IsoD = np.concatenate((np.linspace(1e2,1e3,10), np.linspace(1e3,1e4,10))) #np.linspace(1e3,1e4)
L2norms_T_IsoD, L2norms_dP_dr_IsoD = CT.ConvergenceIsoD(mc_IsoD,me_IsoD,rp_IsoD,pd_IsoD,kd_IsoD,Td_IsoD,lp_IsoD,n_array_IsoD)



#%% L2norms full range 1e2 - 1e5

n_array_IsoD48 = np.concatenate((np.linspace(1e2,1e3,30), np.linspace(1e3,1e4,30), np.linspace(1e4,1e5,30)))
L2norms_T_IsoD48, L2norms_dP_dr_upper_IsoD48 = CT.ConvergenceIsoD(mc_IsoD,me_IsoD,rp_IsoD,pd_IsoD,kd_IsoD,Td_IsoD,lp_IsoD,n_array_IsoD48)



#%% Shows the increase is in fact numerical noise as L2norm per cell average around 1e-7 so for 1e7 cells L2norm is 1e0

n_IsoD2 = 1e7
IsoD2 = IDM.Iso_Density_Giant_Planet(mc_IsoD,me_IsoD,rp_IsoD,pd_IsoD,kd_IsoD,Td_IsoD,lp_IsoD,n_IsoD2)
IsoD2_Sim_Data = IsoD2.Planet_Formation()

L2norms_T2 = IsoD2.L2norm()[5]#, L2norms_dP_dr2 = IsoT2.L2norm()[3:]
avr_L2norm_T2, max_L2norm_T2 = np.mean(L2norms_T2), max(L2norms_T2)
#avr_L2norm_dP_dr2, max_L2norm_dP_dr2 = np.mean(L2norms_dP_dr2), max(L2norms_dP_dr2)
print("Density avr. and max. L2norm =", avr_L2norm_T2, "and", max_L2norm_T2)
#print("Press. grad. avr. and max. L2norm =", avr_L2norm_dP_dr2, "and", max_L2norm_dP_dr2)



#%%


#                       Constant Opacity GPF



#%% Values
# mc, me, rp, pd, kd, Td, lp, n, acc_m, p_ana_c_ini
val_IsoO = [10., 1., 43.5589, 1.47e-11, 6e-5, 118, 0.715809160364261e-4, 1e3, 1e-6, [6.14e-11, 9.20]]
#val_IsoO = [12.07, 1.53, 40., 1e-10, 1e-4, 115, 1e-6, 1e3, 1e-6, [5.38e-11, 9.16]]
mc_IsoO, me_IsoO, rp_IsoO, pd_IsoO, kd_IsoO, Td_IsoO, lp_IsoO, n_IsoO, acc_m, p_ana_c_ini_IsoO = val_IsoO



#%% Calculating the Mass Grid

m_grid_IsoO, x = MGS.Mass_Grid_Solver_Iso_Opacity(mc_IsoO,me_IsoO,rp_IsoO,pd_IsoO,kd_IsoO,Td_IsoO,lp_IsoO,n_IsoO,acc_m,p_ana_c_ini_IsoO)
#m_grid_IsoOM = MGS.Mass_Grid_Solver_Iso_Opacity(mc_IsoO,me_IsoO,rp_IsoO,pd_IsoO,kd_IsoO,Td_IsoO,lp_IsoO,n_IsoO,acc_m,p_ana_c_ini_IsoO)



#%% Running the simulation

IsoO1 = IOM.Iso_Opacity_Giant_Planet(mc_IsoO,me_IsoO,rp_IsoO,pd_IsoO,kd_IsoO,Td_IsoO,lp_IsoO,n_IsoO,grid_given=True,m_grid=m_grid_IsoO)
#IsoOM = IOM.Iso_Opacity_Giant_Planet(mc_IsoO,me_IsoO,rp_IsoO,pd_IsoO,kd_IsoO,Td_IsoO,lp_IsoO,n_IsoO,grid_given=True,m_grid=m_grid_IsoOM)



#%% Running simulation

IsoO1_Sim_Data = IsoO1.Planet_Formation("Both")
#IsoOM_Sim_Data = IsoOM.Planet_Formation("Both")



#%% Plotting Output Data

params = [7, 2, 20, 20, 15, 20]

Plotter.Plot([IsoO1], "m", "r", Params=params)
Plotter.Plot([IsoO1], "p", "r", Params=params)
Plotter.Plot([IsoO1], "k", "r", Params=params) 
Plotter.Plot([IsoO1], "T", "r", Params=params)
Plotter.Plot([IsoO1], "P", "r", Params=params)
Plotter.Plot([IsoO1], "s", "r", Params=params)

#Plotter.Plot([IsoOM], "p", "r", Comparison_Data=Movshovitz_Figure_3)
#Plotter.Plot([IsoOM], "k", "r", Comparison_Data=Movshovitz_Figure_3)
#Plotter.Plot([IsoOM], "T", "r", Comparison_Data=Movshovitz_Figure_3)



#%% Viva Plot
import Plotter
params = [7, 2, 20, 20, 15, 20]

Plotter.Plot([IsoO1], "p+T", "r", Params=params)



#%% Convergence of solution

n_array_IsoO = np.geomspace(1e2,1e5,4)

CT.ConvergenceIsoO(mc_IsoO,me_IsoO,rp_IsoO,pd_IsoO,kd_IsoO,Td_IsoO,lp_IsoO,n_array_IsoO,acc_m,p_ana_c_ini_IsoO)





#%%


#                       Movshovitz Opacity Forced GPF



#%% Values are just test values to show the code works (not necessarily realisitic)
# mc, me, rp, pd, kd, Td, lp, n, acc_m, p_ana_c_ini
val_MovOF3 = [12.07, 1.53, 40., 1e-10, 1., 115, 1e-6, 1e4, 1e-6, [3.28e-11, 1.59e1]]   #[3.28e-11, 1.59e1]     [5.13e-11, 1.08e1]
mc_MovOF3, me_MovOF3, rp_MovOF3, pd_MovOF3, kd_MovOF3, Td_MovOF3, lp_MovOF3, n_MovOF3, acc_m, p_ana_c_ini_MovOF3 = val_MovOF3

#val_MovOF4 = [8., 8., 60., 1e-10, 1e-4, 115, 2e-6, 1e4, 1e-6, [5.13e-11, 1.08e1]]
#mc_MovOF4, me_MovOF4, rp_MovOF4, pd_MovOF4, kd_MovOF4, Td_MovOF4, lp_MovOF4, n_MovOF4, acc_m, p_ana_c_ini_MovOF4 = val_MovOF4



#%% Calculating the Mass Grid

m_grid_MovOF3 = MGS.Mass_Grid_Solver_Movshovitz_Opacity_Forced(mc_MovOF3,me_MovOF3,rp_MovOF3,pd_MovOF3,kd_MovOF3,Td_MovOF3,lp_MovOF3,n_MovOF3,acc_m,p_ana_c_ini_MovOF3,Movshovitz_Figure_3_Opacity_Log)
#m_grid_MovOF4 = MGS.Mass_Grid_Solver_Movshovitz_Opacity_Forced(mc_MovOF4,me_MovOF4,rp_MovOF4,pd_MovOF4,kd_MovOF4,Td_MovOF4,lp_MovOF4,n_MovOF4,acc_m,p_ana_c_ini_MovOF4,Movshovitz_Figure_4_Opacity_Log)



#%% Initialise Giant Planet

MovOF3 = GPFMOF.Movshovitz_Opacity_Forced_Giant_Planet(mc_MovOF3,me_MovOF3,rp_MovOF3,pd_MovOF3,kd_MovOF3,Td_MovOF3,lp_MovOF3,n_MovOF3,p_ana_c=[3.28e-11, 1.59e1])#grid_given=True,m_grid=m_grid_MovOF3)
#MovOF4 = GPFMOF.Movshovitz_Opacity_Forced_Giant_Planet(mc_MovOF4,me_MovOF4,rp_MovOF4,pd_MovOF4,kd_MovOF4,Td_MovOF4,lp_MovOF4,n_MovOF4,p_ana_c=[5.13e-11, 1.08e1])#grid_given=True,m_grid=m_grid_MovOF4)



#%% Running Simulation

MovOF3_Sim_Data = MovOF3.Planet_Formation(Movshovitz_Figure_3_Opacity_Log)
#MovOF4_Sim_Data = MovOF4.Planet_Formation(Movshovitz_Figure_4_Opacity_Log)



#%% Plotting Output Data

MovOF3.Plot("p", "r", Comparison_Data=Movshovitz_Figure_3)
MovOF3.Plot("k", "r", Comparison_Data=Movshovitz_Figure_3)
MovOF3.Plot("T", "r", Comparison_Data=Movshovitz_Figure_3)

#MovOF4.Plot("p", "r", Comparison_Data=Movshovitz_Figure_4)
#MovOF4.Plot("k", "r", Comparison_Data=Movshovitz_Figure_4)
#MovOF4.Plot("T", "r", Comparison_Data=Movshovitz_Figure_4)





#%%


#                       Movshovitz Opacity GPF



#%% Values
# mc, me, rp, pd, kd, Td, lp, n, acc_m, p_ana_c_ini
val_MovO3 = [12.07, 1.53, 40., 4e-11, 1., 115, 1e-6, 1e5, 1e-6, [1, 1]] #p=4e-11
mc_MovO3, me_MovO3, rp_MovO3, pd_MovO3, kd_MovO3, Td_MovO3, lp_MovO3, n_MovO3, acc_m, p_ana_c_ini_MovO3 = val_MovO3

#val_MovO4 = [8., 8., 60., 1e-10, 1e-4, 115, 2e-6, 1e4, 1e-6, [5.13e-11, 1.08e1]]
#mc_MovO4, me_MovO4, rp_MovO4, pd_MovO4, kd_MovO4, Td_MovO4, lp_MovO4, n_MovO4, acc_m, p_ana_c_ini_MovO4 = val_MovO4

"""
#%% Calculating the Mass Grid

m_grid_MovO3 = MGS.Mass_Grid_Solver_Movshovitz_Opacity(mc_MovO3,me_MovO3,rp_MovO3,pd_MovO3,kd_MovO3,Td_MovO3,lp_MovO3,n_MovO3,acc_m,p_ana_c_ini_MovO3,Movshovitz_Figure_3_Opacity_Log)
#m_grid_MovO4 = MGS.Mass_Grid_Solver_Movshovitz_Opacity(mc_MovO4,me_MovO4,rp_MovO4,pd_MovO4,kd_MovO4,Td_MovO4,lp_MovO4,n_MovO4,acc_m,p_ana_c_ini_MovO4,Movshovitz_Figure_4_Opacity_Log)

"""

#%% Initialise Giant Planet

MovO3 = GPFMO.Movshovitz_Opacity_Giant_Planet(mc_MovO3,me_MovO3,rp_MovO3,pd_MovO3,kd_MovO3,Td_MovO3,lp_MovO3,n_MovO3) #(0.7028)
#MovO4 = GPFMO.Movshovitz_Opacity_Giant_Planet(16., 16., 60., 1e-10, 1.70, 115, 3e-6, 1e5, p_ana_c=[7.97989053e-11, 7.40708628]) #(0.7465)



#%% Running Simulation

MovO3_Sim_Data = MovO3.Planet_Formation(Movshovitz_Figure_3_Opacity_Log)
#MovO4_Sim_Data = MovO4.Planet_Formation(Movshovitz_Figure_4_Opacity_Log)



#%% Plotting Output Data

params = [7, 2, 20, 20, 15, 20]

Plotter.Plot([MovO3], "p", "r", Comparison_Data=Movshovitz_Figure_3, Params=params, Movshovitz=True)
Plotter.Plot([MovO3], "k", "r", Comparison_Data=Movshovitz_Figure_3, Params=params, Movshovitz=True)
Plotter.Plot([MovO3], "T", "r", Comparison_Data=Movshovitz_Figure_3, Params=params, Movshovitz=True)
#Plotter.Plot([MovO3], "p+T", "r", Comparison_Data=Movshovitz_Figure_3, Params=params, Movshovitz=True)

#Plotter.Plot([MovO4], "p", "r", Comparison_Data=Movshovitz_Figure_4)
#Plotter.Plot([MovO4], "k", "r", Comparison_Data=Movshovitz_Figure_4)
#Plotter.Plot([MovO4], "T", "r", Comparison_Data=Movshovitz_Figure_4)



#%%% Extra Plots
params = [7, 2, 20, 20, 15, 20]

Plotter.Plot([MovO3], "m", "r", Comparison_Data=Movshovitz_Figure_3, Params=params, Movshovitz=True)
Plotter.Plot([MovO3], "P", "r", Comparison_Data=Movshovitz_Figure_3, Params=params, Movshovitz=True)
Plotter.Plot([MovO3], "s", "r", Comparison_Data=Movshovitz_Figure_3, Params=params, Movshovitz=True)
Plotter.Plot([MovO3], "m", "r", Comparison_Data=Movshovitz_Figure_3, Params=params, Movshovitz=True)



#%% Convergence of solution

n_array_MovO = np.linspace(1e5,1e6,10)

CT.ConvergenceMovO(mc_MovO3,me_MovO3,rp_MovO3,pd_MovO3,kd_MovO3,Td_MovO3,lp_MovO3,n_array_MovO,Movshovitz_Figure_3_Opacity_Log)





#%%


#                       Freedman Opacity GPF

import GPF_Freedman_Opacity as GPFFO

#%% Values
# mc, me, rp, pd, kd, Td, lp, met, acc_k, n, acc_m, p_ana_c_ini
val_FreO = [10., 1., 36.3776, 1e-10, 1e-4, 118, 1e-10, 0., 1e-10, 1e3, 1e-6, [4.69729433e-11, 1.05263989e+01]]

#val_FreO = [10., 1., 43.5589, 1.47e-11, 1, 118, 1e-10, 0., 1e-10, 1e3, 1e-6, [5.58e-11, 9.02]]
#val_FreO = [10., 1., 43.5589, 1.47e-11, 1, 118, 1e-10, 0., 1e-10, 1e3, 1e-6, [4.60123367e+10, 2.71881003e+10]] #36.3776, 38.9459
#val_FreO = [12.07, 1.53, 40., 1e-10, 1e-4, 115, 1e-6, 0., 1e-10, 1e3, 1e-6, [5.58e-11, 9.02]]
mc_FreO, me_FreO, rp_FreO, pd_FreO, kd_FreO, Td_FreO, lp_FreO, met_FreO, acc_k, n_FreO, acc_m, p_ana_c_ini_FreO = val_FreO



#%% Calculating the Mass Grid

m_grid_fre, x = MGS.Mass_Grid_Solver_Freedman_Opacity(mc_FreO,me_FreO,rp_FreO,pd_FreO,kd_FreO,Td_FreO,lp_FreO,met_FreO,acc_k,n_FreO,acc_m,p_ana_c_ini_FreO)



#%% Initialise Giant Planet

#FreO1 = GPFFO.Freedman_Opacity_Giant_Planet(mc_FreO,me_FreO,rp_FreO,1.47e-11,kd_FreO,Td_FreO,4.4764e-06,met_FreO,acc_k,n_FreO,grid_given=True,m_grid=m_grid_fre)
FreO1 = GPFFO.Freedman_Opacity_Giant_Planet(mc_FreO,me_FreO,rp_FreO,pd_FreO,kd_FreO,Td_FreO,lp_FreO,met_FreO,acc_k,n_FreO,grid_given=True,m_grid=m_grid_fre)
#FreOM = GPFFO.Freedman_Opacity_Giant_Planet(mc_FreO,me_FreO,rp_FreO,pd_FreO,kd_FreO,Td_FreO,lp_FreO,met_FreO,acc_k,n_FreO,grid_given=True,m_grid=m_grid_fre)



#%% Running Simulation

FreO1_Sim_Data = FreO1.Planet_Formation()
#FreOM_Sim_Data = FreOM.Planet_Formation()



#%% Plotting Output Data

params = [7, 2, 20, 20, 15, 20]

Plotter.Plot([FreO1], "m", "r", Params=params)
Plotter.Plot([FreO1], "p", "r", Params=params)
Plotter.Plot([FreO1], "k", "r", Params=params) 
Plotter.Plot([FreO1], "T", "r", Params=params)
Plotter.Plot([FreO1], "P", "r", Params=params)
Plotter.Plot([FreO1], "s", "r", Params=params)

#Plotter.Plot([FreOM], "p", "r", Comparison_Data=Movshovitz_Figure_3)
#Plotter.Plot([FreOM], "k", "r", Comparison_Data=Movshovitz_Figure_3)
#Plotter.Plot([FreOM], "T", "r", Comparison_Data=Movshovitz_Figure_3)



#%% Viva Plot
import Plotter
params = [7, 2, 20, 20, 15, 20]

Plotter.Plot([FreO1], "k+p+T", "r", Params=params)



#%% Convergence of solution

n_array_FreO = np.linspace(1e2,1e3,10)

CT.ConvergenceFreO(mc_FreO,me_FreO,rp_FreO,pd_FreO,kd_FreO,Td_FreO,lp_FreO,met_FreO,acc_k,n_array_FreO,acc_m,p_ana_c_ini_FreO)





#%%


#                       Combined Opacity GPF



#%% Uniform logarithmic spacing for me (envelope mass)
M_E = 5.972e27
env_mass_min = 1e20
env_mass_max = 1. * M_E
num_p = 1e3
m1 = np.geomspace(env_mass_min, 0.9*env_mass_max, num=int(0.3*num_p))
m2 = np.geomspace((0.9*env_mass_max), env_mass_max, num=int((0.7*num_p)+1))
m2 = m2[1:]
m = np.concatenate((m1,m2))

for i in range(1):
    m = np.insert(m, 0,      0.)
    m = np.insert(m, len(m), 0.)
m_grid_com = m



#%% Values
# mc, me, rp, pd, kd, Td, lp, ad, d_to_g, acc_k, n, acc_m, p_ana_c_ini
val_ComO = [10., 1., 30., 1e-10, 1e-4, 115, 1e-6, 1e-2, 0.01, 1e-10, 1e3, 1e-6, [1, 1]]
mc_ComO, me_ComO, rp_ComO, pd_ComO, kd_ComO, Td_ComO, lp_ComO, ad_ComO, d_to_g_ComO, acc_k, n_ComO, acc_m, p_ana_c_ini_ComO = val_ComO



#%% Calculating the Mass Grid
import Mass_Grid_Solver_Combined_Opacity as MGSCO
# mc, me, rp, pd, kd, Td, lp, ad, d_to_g, acc_k, n, acc_m, p_ana_c_ini
m_grid_com = MGSCO.Mass_Grid_Solver(mc_ComO,me_ComO,rp_ComO,pd_ComO,kd_ComO,Td_ComO,lp_ComO,ad_ComO,d_to_g_ComO,acc_k,n_ComO,acc_m,p_ana_c_ini_ComO)



#%% Initialise Giant Planet
# mc, me, rp, pd, kd, Td, lp, ad, d_to_g, acc_k, n
#ComO1 = GPFCO.Combined_Opacity_Giant_Planet(10., 1., 30., 1e-10, 1e-4, 115, 1e-6, ad, d_to_g, 1e-10, 1e4, p_ana_c=[5.08e-6, -1.42e6])# grid_given=True, m_grid=m_grid_com)
ComO1 = GPFCO.Combined_Opacity_Giant_Planet(mc_ComO,me_ComO,rp_ComO,pd_ComO,kd_ComO,Td_ComO,lp_ComO,ad_ComO,d_to_g_ComO,acc_k,n_ComO,grid_given=True,m_grid=m_grid_com)



#%% Running Simulation

ComO1_Sim_Data = ComO1.Planet_Formation()



#%% Plotting Output Data

Plotter.Plot([ComO1], "p", "r", Comparison_Data=Movshovitz_Figure_3)
Plotter.Plot([ComO1], "k", "r", Comparison_Data=Movshovitz_Figure_3)
Plotter.Plot([ComO1], "T", "r", Comparison_Data=Movshovitz_Figure_3)



#%% Convergence of solution

n_array_ComO = np.linspace(1e3,1e4,10)

CT.ConvergenceComO(mc_ComO,me_ComO,rp_ComO,pd_ComO,kd_ComO,Td_ComO,lp_ComO,ad_ComO,d_to_g_ComO,acc_k,n_array_ComO,acc_m,p_ana_c_ini_ComO)





#%%


#                       Maximum Mass Accretion Testing


import Maximum_Accretion_Rate_Calculator as MARC
#%% Values
# Sigma_outer, Ms, Mp, Tp, Rp, r_in, r_out, alpha, num_sig
val_Sigma = [143.365/np.sqrt(3), 1., 11.*5.972e27, 118., 5.2, 2., 3., 1e-5, 1e4]
Val_outer_Sigma, Ms_Sigma, Mp_Sigma, Tp_Sigma, Rp_Sigma, r_in_Sigma, r_out_Sigma, alpha_Sigma, num_sig_Sigma = val_Sigma



#%% Testing 0 mass planet

# Sigma_outer, Ms, Mp, Tp, Rp, r_in, r_out, alpha, num_sig
Mdot_max, Sigma, r = MARC.Max_Acc_Rate(Val_outer_Sigma,Ms_Sigma,0.,Tp_Sigma,Rp_Sigma,r_in_Sigma,r_out_Sigma,alpha_Sigma,num_sig_Sigma,Plot=True)

# No gap forms in disc as expected and -> perfect 1/sqrt(r) dependence in Sigma



#%% Testing Crida Sigma eqn 14 Sigma solution - comparing to Figure 11
import Maximum_Accretion_Rate_Calculator as MARC
# Crida Comparison -> Hard sets: Sigma_outer = 1/np.sqrt(3), nu = 3e14, H/r = 0.05, Mp = 1e-3 * Ms
# Sigma_outer, Ms, Mp, Tp, Rp, r_in, r_out, alpha, num_sig
Sigma, R_sig = MARC.Max_Acc_Rate(1, 1., 1, 1, 5.2, 2., 3., 1, 1e5, Plot=True, Crida_Comparison=True)



#%% Different alpha values effect

import matplotlib.pyplot as plt
import scipy.constants as spc


AU = spc.au * 1e2
alphas = [1e-2, 1e-5]

_, Sigma0, r0 = MARC.Max_Acc_Rate(Val_outer_Sigma,Ms_Sigma,11.*5.972e27,Tp_Sigma,Rp_Sigma,r_in_Sigma,r_out_Sigma,alphas[0],num_sig_Sigma)
_, Sigma1, r1 = MARC.Max_Acc_Rate(Val_outer_Sigma,Ms_Sigma,11.*5.972e27,Tp_Sigma,Rp_Sigma,r_in_Sigma,r_out_Sigma,alphas[1],num_sig_Sigma)
_, Sigma2, r2 = MARC.Max_Acc_Rate(Val_outer_Sigma,Ms_Sigma,300.*5.972e27,Tp_Sigma,Rp_Sigma,r_in_Sigma,r_out_Sigma,alphas[0],num_sig_Sigma)
_, Sigma3, r3 = MARC.Max_Acc_Rate(Val_outer_Sigma,Ms_Sigma,300.*5.972e27,Tp_Sigma,Rp_Sigma,1,1,alphas[1],num_sig_Sigma,L09=True)

darkness = 0.5

Figsize, lw, Ticksize, Labelsize, Legendfontsize = [7, 2, 20, 20, 15]

plt.rc("figure", figsize=(Figsize,Figsize))
plt.rc("xtick", labelsize=Ticksize)
plt.rc("ytick", labelsize=Ticksize)
plt.rc("axes", labelsize=Labelsize)
plt.rc("legend", fontsize=Legendfontsize)

plt.figure()
plt.plot(r0/AU, Sigma0, linestyle="-", color="blue", label=r"$\alpha$=10$^{-2}$, 0.1 $M_{\oplus}$")
plt.plot(r2/AU, Sigma2, linestyle="--", color="blue", label=r"$\alpha$=10$^{-2}$, 300 $M_{\oplus}$")
plt.plot(r1/AU, Sigma1, linestyle="-", color="orange", label=r"$\alpha$=10$^{-5}$, 0.1 $M_{\oplus}$")
plt.xlabel("Radius / AU")
plt.ylabel("$\Sigma$ / g$\cdot$cm$^{-2}$")
plt.legend()
ax = plt.gca()
ax.set_box_aspect(1)
plt.tight_layout()
plt.show()

plt.figure()
plt.plot(r2/AU, Sigma2, linestyle="-", color="blue", label=r"$\alpha$=10$^{-2}$, 300$M_{\oplus}$")
plt.axvspan(r2[int(num_sig_Sigma/5)]/AU, r2[int(2*num_sig_Sigma/5)-1]/AU, alpha=darkness, facecolor="red", label="T12")
plt.axvspan(r2[int(3*num_sig_Sigma/5)-1]/AU, r2[int(4*num_sig_Sigma/5)-1]/AU, alpha=darkness, facecolor="red")
plt.axvspan(r3[int(num_sig_Sigma/5)]/AU, r3[int(2*num_sig_Sigma/5)-1]/AU, alpha=darkness, facecolor="orange", label="L09")
plt.axvspan(r3[int(3*num_sig_Sigma/5)-1]/AU, r3[int(4*num_sig_Sigma/5)-1]/AU, alpha=darkness, facecolor="orange")
plt.xlabel("Radius / AU")
plt.ylabel("$\Sigma$ / g$\cdot$cm$^{-2}$")
plt.legend()
ax = plt.gca()
ax.set_box_aspect(1)
plt.tight_layout()
plt.show()



#%% Convergence of Sigma solution
import Convergence_Tester as CT
n_array_Sigma = np.geomspace(1e2,1e5,4)

CT.ConvergenceSigma(Val_outer_Sigma,Ms_Sigma,Mp_Sigma,Tp_Sigma,Rp_Sigma,r_in_Sigma,r_out_Sigma,alpha_Sigma,n_array_Sigma)




#%% Convergence of Mdot solution with analytic solution

n_array_Sigma = np.geomspace(1e2,1e5,4)

CT.ConvergenceMdotMax(Val_outer_Sigma,Ms_Sigma,Mp_Sigma,Tp_Sigma,Rp_Sigma,r_in_Sigma,r_out_Sigma,alpha_Sigma,n_array_Sigma)

# Can be confident that the Mdot_max equation is suitably implemented





#%%


#                       Luminosity Finder IsoO Testing



#%% Values
# mc, me, rp, pd, kd, Td, lp, acc_l, n, acc_m, p_ana_c_ini
val_IsoO = [10., 1., 43.5589, 1.47e-11, 6e-5, 118, 0.72e-4, 1e-10, 1e3, 1e-6, [6.14e-11, 9.20]] #38.9459
mc_IsoO, me_IsoO, rp_IsoO, pd_IsoO, kd_IsoO, Td_IsoO, lp_IsoO, acc_l, n_IsoO, acc_m, p_ana_c_ini_IsoO = val_IsoO

import L_Finder_IsoO as LFIO

#%% Calculating the initial Mass Grid

m_grid_IsoO, x = MGS.Mass_Grid_Solver_Iso_Opacity(mc_IsoO,me_IsoO,rp_IsoO,pd_IsoO,kd_IsoO,Td_IsoO,lp_IsoO,n_IsoO,acc_m,p_ana_c_ini_IsoO)



#%% Finding L for given inner radius

#L1 = LFIO.L_Finder(mc_IsoO,me_IsoO,rp_IsoO,pd_IsoO,kd_IsoO,Td_IsoO,lp_IsoO,acc_l,n_IsoO,m_grid_IsoO)[0]
L1 = LFIO.L_Finder(mc_IsoO,me_IsoO,rp_IsoO,pd_IsoO,kd_IsoO,Td_IsoO,lp_IsoO,acc_l,n_IsoO,m_grid_IsoO)[0]
print(L1)



#%% Convergence Test on L Finding

CT.ConvergenceLumIsoO(mc_IsoO,me_IsoO,rp_IsoO,pd_IsoO,kd_IsoO,Td_IsoO,lp_IsoO,acc_l,n_IsoO,m_grid_IsoO)





#%%


#                           Evolver IsoO

import Evolver_IsoO as EIO
import Plotter
#%% Values
# mc, me, rp, pd, kd, Td, lp, acc_l, n_i, acc_m, p_ana_c_ini, t_tot, n_t
val_IsoO = [10., 1., 43.5589, 1.47e-11, 6e-5, 118, 0.72e-4, 1e-10, 1e4, 1e-6, [6.14e-11, 9.20], 2e6, 1000] #[7.58e-11, 8.65], [7.5614e-11, 8.6715], 38.9459
mc_IsoO, me_IsoO, rp_IsoO, pd_IsoO, kd_IsoO, Td_IsoO, lp_IsoO, acc_l, n_i_IsoO, acc_m, p_ana_c_ini_IsoO, t_tot_IsoO, n_t_IsoO = val_IsoO
# 7.2e-11, 700
# Sigma_outer, Ms, Rp, r_in, r_out, alpha, num_sig
val_Sigma = [143.365/np.sqrt(3), 1., 5.2, 2, 3, 1e-5, 1e4]
Val_outer_IsoO, Ms_IsoO, Rp_IsoO, r_in_IsoO, r_out_IsoO, alpha_IsoO, num_sig_IsoO = val_Sigma



#%% Calculating the initial Mass Grid

m_grid_IsoO, _ = MGS.Mass_Grid_Solver_Iso_Opacity(mc_IsoO,me_IsoO,rp_IsoO,pd_IsoO,kd_IsoO,Td_IsoO,lp_IsoO,n_i_IsoO,acc_m,p_ana_c_ini_IsoO)



#%%

start_time = time.time()

#EIO1 = EIO.EvolverIsoO(mc_IsoO,me_IsoO,rp_IsoO,pd_IsoO,kd_IsoO,Td_IsoO,lp_IsoO,acc_l,n_i_IsoO,acc_m,p_ana_c_ini_IsoO,m_grid_IsoO,t_tot_IsoO,n_t_IsoO,Val_outer_IsoO,Ms_IsoO,Rp_IsoO,r_in_IsoO,r_out_IsoO,alpha_IsoO,num_sig_IsoO,L09=True)
EIO2 = EIO.EvolverIsoO(mc_IsoO,me_IsoO,rp_IsoO,pd_IsoO,kd_IsoO,Td_IsoO,lp_IsoO,acc_l,n_i_IsoO,acc_m,p_ana_c_ini_IsoO,m_grid_IsoO,t_tot_IsoO,n_t_IsoO,Val_outer_IsoO,Ms_IsoO,Rp_IsoO,r_in_IsoO,r_out_IsoO,alpha_IsoO,num_sig_IsoO,L09=True)

end_time = time.time()
print("Runtime =", end_time - start_time, "seconds")



#%% Calc for plots
params = [7, 2, 20, 20, 15, 20]
import Plotter
#Calc1 = Plotter.Plot([EIO1], _, "calc", Params=params)
Calc2 = Plotter.Plot([EIO2], _, "calc", Params=params)



#%% Plots
params = [7, 2, 20, 20, 15, 20]
import Plotter
#Plotter.Plot([EIO1], "M", "t", Params=params, calc=Calc1)
#Plotter.Plot([EIO1], "L", "t", Params=params, calc=Calc1)
#Plotter.Plot([EIO1], "L_t", "t", Params=params, calc=Calc1)
#Plotter.Plot([EIO1], "T_t,r_t", "t", Params=params, calc=Calc1)
Plotter.Plot([EIO1], "T_t4,r_t2", "t", Params=params, calc=Calc1)
#Plotter.Plot([EIO1], "T_t,p_t", "r_t", Params=params, Evo=True)
#Plotter.Plot([EIO1], "p_t comp", "r_t", Params=params, Evo=True)
#Plotter.Plot([EIO1], "m_t comp", "r_t", Params=params, Evo=True)
#Plotter.Plot([EIO1], "T_t,p_t diff", "r_t", Params=params, Evo=True)
#Plotter.Plot([EIO1], "T_t", "r_t", Params=params, Evo=True)
Plotter.Plot([EIO1], "L+M", "t", Params=params, calc=Calc1)
#Plotter.Plot([EIO1], "Mdot", "t", Params=params, calc=Calc1)



#%% Plots
params = [7, 2, 20, 20, 15, 20]
import Plotter
#Plotter.Plot([EIO2], "M", "t", Params=params, calc=Calc2)
#Plotter.Plot([EIO2], "L", "t", Params=params, calc=Calc2)
#Plotter.Plot([EIO2], "L_t", "t", Params=params, calc=Calc2)
#Plotter.Plot([EIO2], "T_t,r_t", "t", Params=params, calc=Calc2)
#Plotter.Plot([EIO2], "T_t4,r_t2", "t", Params=params, calc=Calc2)
#Plotter.Plot([EIO2], "T_t,p_t", "r_t", Params=params, Evo=True)
#Plotter.Plot([EIO2], "p_t comp", "r_t", Params=params, Evo=True)
Plotter.Plot([EIO2], "m_t comp", "r_t", Params=params, Evo=True)
#Plotter.Plot([EIO2], "T_t,p_t diff", "r_t", Params=params, Evo=True)
#Plotter.Plot([EIO2], "T_t", "r_t", Params=params, Evo=True)
#Plotter.Plot([EIO2], "L+M", "t", Params=params, calc=Calc2)
#Plotter.Plot([EIO2], "Mdot", "t", Params=params, calc=Calc2)



#%% Different constant opacity evolutions

kd_IsoOA, kd_IsoOB, kd_IsoOC = 6e-5, 1e-5, 0.5e-5
lp_IsoOA, lp_IsoOB, lp_IsoOC = 0.70e-4, 4.14e-4, 8.29e-4

#%%
m_grid_IsoOA, xA = MGS.Mass_Grid_Solver_Iso_Opacity(mc_IsoO,me_IsoO,rp_IsoO,pd_IsoO,kd_IsoOA,Td_IsoO,lp_IsoOA,n_i_IsoO,acc_m,p_ana_c_ini_IsoO)
EIOA = EIO.EvolverIsoO(mc_IsoO,me_IsoO,rp_IsoO,pd_IsoO,kd_IsoOA,Td_IsoO,lp_IsoOA,acc_l,n_i_IsoO,acc_m,p_ana_c_ini_IsoO,m_grid_IsoOA,t_tot_IsoO,n_t_IsoO,Val_outer_IsoO,Ms_IsoO,Rp_IsoO,r_in_IsoO,r_out_IsoO,alpha_IsoO,num_sig_IsoO)
#%%
m_grid_IsoOB, xB = MGS.Mass_Grid_Solver_Iso_Opacity(mc_IsoO,me_IsoO,rp_IsoO,pd_IsoO,kd_IsoOB,Td_IsoO,lp_IsoOB,n_i_IsoO,acc_m,p_ana_c_ini_IsoO)
EIOB = EIO.EvolverIsoO(mc_IsoO,me_IsoO,rp_IsoO,pd_IsoO,kd_IsoOB,Td_IsoO,lp_IsoOB,acc_l,n_i_IsoO,acc_m,p_ana_c_ini_IsoO,m_grid_IsoOB,t_tot_IsoO,n_t_IsoO,Val_outer_IsoO,Ms_IsoO,Rp_IsoO,r_in_IsoO,r_out_IsoO,alpha_IsoO,num_sig_IsoO)
#%%
m_grid_IsoOC, xC = MGS.Mass_Grid_Solver_Iso_Opacity(mc_IsoO,me_IsoO,rp_IsoO,pd_IsoO,kd_IsoOC,Td_IsoO,lp_IsoOC,3*n_i_IsoO,acc_m,p_ana_c_ini_IsoO)
EIOC = EIO.EvolverIsoO(mc_IsoO,me_IsoO,rp_IsoO,pd_IsoO,kd_IsoOC,Td_IsoO,lp_IsoOC,acc_l,3*n_i_IsoO,acc_m,p_ana_c_ini_IsoO,m_grid_IsoOC,t_tot_IsoO,n_t_IsoO,Val_outer_IsoO,Ms_IsoO,Rp_IsoO,r_in_IsoO,r_out_IsoO,alpha_IsoO,num_sig_IsoO)



#%% Plotting different constant opacity evolutions

params = [7, 2, 20, 20, 15, 20]
import Plotter

Plotter.Plot([EIOA,EIOB,EIOC], "M", "t_pres", Params=params, Evo=True)
Plotter.Plot([EIOA,EIOB,EIOC], "L", "t_pres", Params=params, Evo=True)



#%% Different alpha values evolutions

alpha_IsoOD, alpha_IsoOE = 1e-5, 1e-2



#%%
m_grid_IsoOD, xD = MGS.Mass_Grid_Solver_Iso_Opacity(mc_IsoO,me_IsoO,rp_IsoO,pd_IsoO,kd_IsoO,Td_IsoO,lp_IsoO,n_i_IsoO,acc_m,p_ana_c_ini_IsoO)
EIOD = EIO.EvolverIsoO(mc_IsoO,me_IsoO,rp_IsoO,pd_IsoO,kd_IsoO,Td_IsoO,lp_IsoO,acc_l,n_i_IsoO,acc_m,p_ana_c_ini_IsoO,m_grid_IsoO,t_tot_IsoO,n_t_IsoO,Val_outer_IsoO,Ms_IsoO,Rp_IsoO,r_in_IsoO,r_out_IsoO,alpha_IsoOD,num_sig_IsoO)
#%%
m_grid_IsoOE, xE = MGS.Mass_Grid_Solver_Iso_Opacity(mc_IsoO,me_IsoO,rp_IsoO,pd_IsoO,kd_IsoO,Td_IsoO,lp_IsoO,n_i_IsoO,acc_m,p_ana_c_ini_IsoO)
EIOE = EIO.EvolverIsoO(mc_IsoO,me_IsoO,rp_IsoO,pd_IsoO,kd_IsoO,Td_IsoO,lp_IsoO,acc_l,n_i_IsoO,acc_m,p_ana_c_ini_IsoO,m_grid_IsoO,t_tot_IsoO,n_t_IsoO,Val_outer_IsoO,Ms_IsoO,Rp_IsoO,r_in_IsoO,r_out_IsoO,alpha_IsoOE,num_sig_IsoO)
#%%
m_grid_IsoOF, xF = MGS.Mass_Grid_Solver_Iso_Opacity(mc_IsoO,me_IsoO,rp_IsoO,pd_IsoO,kd_IsoO,Td_IsoO,lp_IsoO,n_i_IsoO,acc_m,p_ana_c_ini_IsoO)
EIOF = EIO.EvolverIsoO(mc_IsoO,me_IsoO,rp_IsoO,pd_IsoO,kd_IsoO,Td_IsoO,lp_IsoO,acc_l,n_i_IsoO,acc_m,p_ana_c_ini_IsoO,m_grid_IsoO,t_tot_IsoO,n_t_IsoO,Val_outer_IsoO,Ms_IsoO,Rp_IsoO,r_in_IsoO,r_out_IsoO,alpha_IsoOD,num_sig_IsoO,L09=True)



#%% Plotting different alpha values evolutions

params = [7, 2, 20, 20, 12, 20]
import Plotter

Plotter.Plot([EIOD,EIOE,EIOF], "Mdot", "t_pres", Params=params, Evo=True, acc_reg=True)



#%% Different constant opacity evolutions

m_grid_IsoOG, xG = MGS.Mass_Grid_Solver_Iso_Opacity(mc_IsoO,me_IsoO,rp_IsoO,pd_IsoO,kd_IsoOA,Td_IsoO,lp_IsoOA,n_i_IsoO,acc_m,p_ana_c_ini_IsoO)
EIOG = EIO.EvolverIsoO(mc_IsoO,me_IsoO,rp_IsoO,pd_IsoO,kd_IsoOA,Td_IsoO,lp_IsoOA,acc_l,n_i_IsoO,acc_m,p_ana_c_ini_IsoO,m_grid_IsoOG,t_tot_IsoO,n_t_IsoO,Val_outer_IsoO,Ms_IsoO,Rp_IsoO,r_in_IsoO,r_out_IsoO,alpha_IsoOD,num_sig_IsoO,L09=True)
#%%
m_grid_IsoOH, xH = MGS.Mass_Grid_Solver_Iso_Opacity(mc_IsoO,me_IsoO,rp_IsoO,pd_IsoO,kd_IsoOB,Td_IsoO,lp_IsoOB,n_i_IsoO,acc_m,p_ana_c_ini_IsoO)
EIOH = EIO.EvolverIsoO(mc_IsoO,me_IsoO,rp_IsoO,pd_IsoO,kd_IsoOB,Td_IsoO,lp_IsoOB,acc_l,n_i_IsoO,acc_m,p_ana_c_ini_IsoO,m_grid_IsoOH,t_tot_IsoO,n_t_IsoO,Val_outer_IsoO,Ms_IsoO,Rp_IsoO,r_in_IsoO,r_out_IsoO,alpha_IsoOD,num_sig_IsoO,L09=True)
#%%
m_grid_IsoOI, xI = MGS.Mass_Grid_Solver_Iso_Opacity(mc_IsoO,me_IsoO,rp_IsoO,pd_IsoO,kd_IsoOC,Td_IsoO,lp_IsoOC,3*n_i_IsoO,acc_m,p_ana_c_ini_IsoO)
EIOI = EIO.EvolverIsoO(mc_IsoO,me_IsoO,rp_IsoO,pd_IsoO,kd_IsoOC,Td_IsoO,lp_IsoOC,acc_l,3*n_i_IsoO,acc_m,p_ana_c_ini_IsoO,m_grid_IsoOI,t_tot_IsoO,n_t_IsoO,Val_outer_IsoO,Ms_IsoO,Rp_IsoO,r_in_IsoO,r_out_IsoO,alpha_IsoOD,num_sig_IsoO,L09=True)



#%% Plotting different constant opacity evolutions

params = [7, 2, 20, 20, 15, 20]
import Plotter

Plotter.Plot([EIOG,EIOH,EIOI], "Mdot", "t_pres", Params=params, Evo=True)



#%% Plotting mdot and m evolution

params = [7, 2, 20, 20, 15, 20]
import Plotter

Plotter.Plot([EIOI], "Mdot+M", "t", Params=params)





#%%


#                       Luminosity Finder Freedman Testing



#%% Values
# mc, me, rp, pd, kd, Td, lp, acc_l, met, acc_k, n, acc_m, p_ana_c_ini

36.3776, 1.47e-11
val_FreO = [10., 1., 36.3776, 1e-10, 1e-4, 118, 4.5e-6, 1e-10, 0., 1e-10, 1e3, 1e-6, [4.69729433e-11, 1.05263989e+01]]
mc_FreO, me_FreO, rp_FreO, pd_FreO, kd_FreO, Td_FreO, lp_FreO, acc_l, met_FreO, acc_k, n_FreO, acc_m, p_ana_c_ini_FreO = val_FreO



#%% Calculating the Mass Grid

m_grid_FreO, _ = MGS.Mass_Grid_Solver_Freedman_Opacity(mc_FreO,me_FreO,rp_FreO,pd_FreO,kd_FreO,Td_FreO,lp_FreO,met_FreO,acc_k,n_FreO,acc_m,p_ana_c_ini_FreO)



#%% Finding L for given inner radius

L1 = LFFO.L_Finder(mc_FreO,me_FreO,rp_FreO,pd_FreO,kd_FreO,Td_FreO,lp_FreO,acc_l,met_FreO,acc_k,n_FreO, m_grid_FreO)[0]
print(L1)



#%% Calculating the initial mass grid

# mc, me, rp, pd, kd, Td, lp, met, acc_k, n, acc_m, p_ana_c_ini
m_grid = MGS.Mass_Grid_Solver_Freedman_Opacity(10., 1., 30., 1e-10, 1, 118, 1e-6, 0., 1e-10, 1e3, 1e-6, [7.59e-11, 8.68])


#%% Finding L for given inner radius

start_time = time.time()

# mc, me, rp, pd, kd, Td, lp, Z, acc_k, p_ana_c_ini, acc, n
L1, var1, FreO1 = LFFO.L_Finder(10., 1., 36.3776, 1e-10, 1e-4, 118, 6.88e-6, 0., 1e-10, 1e-10, 1e3, m_grid_fre)[:3]

end_time = time.time()

print("Runtime =", end_time - start_time, "seconds")
print(L1)



#%% Convergence Test on L Finding

# mc, me, rp, pd, kd, Td, lp, Z, acc_k, p_ana_c_ini, rc, acc, n
CT.ConvergenceLumFreO(10., 1., 30., 1e-10, 1e-4, 115, 1e-6, 0., 1e-10, 3.6e9, 1e-10, 1e3, m_grid)





#%%


#                           Evolver Freedman



#%% Calculating the Initial Mass Grid

# mc, me, rp, pd, kd, Td, lp, met, acc_k, n, acc_m, p_ana_c_ini
m_grid = MGS.Mass_Grid_Solver_Freedman_Opacity(10., 1., 30., 1e-10, 1e-4, 115, 1e-6, 0., 1e-10, 1e3, 1e-6, [7.59e-11, 8.68])



#%%

start_time = time.time()

# mc, me, rp, pd, kd, Td, L_guess_ini, Z, acc_k, rc, acc, n_i, m_grid, t_tot, n_t, Sigma_in, Ms, Rp, r_in, r_out, alpha, num_sig
L1, M1, t1, varis1 = EFO.Evolver(10., 1., 30., 1e-10, 1e-4, 115, 1e-6, 0., 1e-10, 3.6e9, 1e-10, 1e3, m_grid, 1e6, 10, 800/np.sqrt(3), 1., 5.2, 2., 3., 1e-3, 1e5)

end_time = time.time()
print("Runtime =", end_time - start_time, "seconds")



#%% Plots

Plotter.Plot([[L1, M1, t1, varis1]], "M", "t")
Plotter.Plot([[L1, M1, t1, varis1]], "L", "t")





#%%


#                       Luminosity Finder Combined Testing

ad, d_to_g = 1e-4, 0.01

#%% Calculating the Initial Mass Grid
import Mass_Grid_Solver_Combined_Opacity as MGSCO
# mc, me, rp, pd, kd, Td, lp, ad, d_to_g, acc_k, n, acc_m, p_ana_c_ini
m_grid = MGSCO.Mass_Grid_Solver(12.07, 1.53, 40., 1e-10, 1e-4, 115, 1e-6, ad, d_to_g, 1e-10, 1e3, 1e-6, [5.59e-11, 9.01])



#%% Finding L for given inner radius

start_time = time.time()
import GPF_Combined_Opacity as GPFCO
import L_Finder_ComO as LFCO
# mc, me, rp, pd, kd, Td, lp, ad, d_to_g, acc_k, p_ana_c_ini, rc, acc, n
L1, var1, ComO1 = LFCO.L_Finder(12.07, 1.53, 40., 1e-10, 1e-4, 115, 1e-6, ad, d_to_g, 1e-10, 3.6e9, 1e-10, 1e3, m_grid)[:3]

end_time = time.time()

print("Runtime =", end_time - start_time, "seconds")
print(L1)



#%% Convergence Test on L Finding

# mc, me, rp, pd, kd, Td, lp, ad, d_to_g, acc_k, p_ana_c_ini, rc, acc, n
CT.ConvergenceLumComO(12.07, 1.53, 40., 1e-10, 1e-4, 115, 1e-6, ad, d_to_g, 1e-10, 3.6e9, 1e-10, 1e3, m_grid)





#%%


#                           Evolver Combined



#%% Calculating the Initial Mass Grid
import Mass_Grid_Solver_Combined_Opacity as MGSCO
# mc, me, rp, pd, kd, Td, lp, ad, d_to_g, acc_k, n, acc_m, p_ana_c_ini
m_grid = MGSCO.Mass_Grid_Solver(12.07, 1.53, 40., 1e-10, 1e-4, 115, 1e-6, ad, d_to_g, 1e-10, 1e3, 1e-6, [5.59e-11, 9.01])



#%%

start_time = time.time()
import GPF_Combined_Opacity as GPFCO
import L_Finder_ComO as LFCO
import Evolver_ComO as ECO
# mc, me, rp, pd, kd, Td, L_guess_ini, ad, d_to_g, acc_k, rc, acc, n_i, m_grid, t_tot, n_t, Sigma_in, Ms, Rp, r_in, r_out, alpha, num_sig
L1, M1, t1, varis1 = ECO.Evolver(12.07, 1.53, 40., 1e-10, 1e-4, 115, 1e-6, ad, d_to_g, 1e-10, 3.6e9, 1e-10, 1e3, m_grid, 1e5, 10, 10, 1., 5.2, 2, 3, 1e-3, 1e3)

end_time = time.time()
print("Runtime =", end_time - start_time, "seconds")



#%% Plots

Plotter.Plot([[L1, M1, t1, varis1]], "M", "t")
Plotter.Plot([[L1, M1, t1, varis1]], "L", "t")





#%%


#                               Poster Plots



#%% Values
# mc, me, rp, pd, kd, Td, lp, n, acc_m, p_ana_c_ini
val_MovO3 = [12.07, 1.53, 40., 1e-10, 1e-4, 115, 1e-6, 1e5, 1e-6]
mc_MovO3, me_MovO3, rp_MovO3, pd_MovO3, kd_MovO3, Td_MovO3, lp_MovO3, n_MovO3, acc_m = val_MovO3



#%% Structure Simulations

# mc, me, rp, pd, kd, Td, lp, n
IsoO1 = IOM.Iso_Opacity_Giant_Planet(mc_MovO3,me_MovO3,rp_MovO3,pd_MovO3,kd_MovO3,Td_MovO3,lp_MovO3,n_MovO3,p_ana_c=[7.25e-11, 8.46])
IsoO1_Sim_Data = IsoO1.Planet_Formation("Both")

# mc, me, rp, pd, kd, Td, lp, n
MovO1 = GPFMO.Movshovitz_Opacity_Giant_Planet(mc_MovO3,me_MovO3,rp_MovO3,pd_MovO3,kd_MovO3,Td_MovO3,lp_MovO3,n_MovO3)
MovO1_Sim_Data = MovO1.Planet_Formation(Movshovitz_Figure_3_Opacity_Log)

# mc, me, rp, pd, kd, Td, lp, Z, acc_k, p_ana_c, n
FreO1 = GPFFO.Freedman_Opacity_Giant_Planet(mc_MovO3,me_MovO3,rp_MovO3,pd_MovO3,kd_MovO3,Td_MovO3,12.467*lp_MovO3,0.,1e-10,n_MovO3,p_ana_c=[5.58e-11, 9.02])
FreO1_Sim_Data = FreO1.Planet_Formation()

# mc, me, rp, pd, kd, Td, lp, Z, acc_k, p_ana_c, n
FreO2 = GPFFO.Freedman_Opacity_Giant_Planet(mc_MovO3,me_MovO3,rp_MovO3,pd_MovO3,kd_MovO3,Td_MovO3,lp_MovO3,0.,1e-10,n_MovO3,p_ana_c=[7.59e-11, 8.68])
FreO2_Sim_Data = FreO2.Planet_Formation()


IsoO1._boundary_i.append([-0,IsoO1._conv_or_rad[-1]])
IsoO1._boundary_i.insert(0, [1,IsoO1._conv_or_rad[0]])
MovO1._boundary_i.append([-0,MovO1._conv_or_rad[-1]])
MovO1._boundary_i.insert(0, [1,MovO1._conv_or_rad[0]])
FreO1._boundary_i.append([-0,FreO1._conv_or_rad[-1]])
FreO1._boundary_i.insert(0, [1,FreO1._conv_or_rad[0]])
FreO2._boundary_i.append([-0,FreO2._conv_or_rad[-1]])
FreO2._boundary_i.insert(0, [1,FreO2._conv_or_rad[0]])


# mc, me, rp, pd, kd, Td, lp, n
#IsoO2_Sim_Data = IsoO2.Planet_Formation("Both")
#IsoO2 = IOM.Iso_Opacity_Giant_Planet(12.07, 1.53, 40., 1e-10, 1e-3, 115, 1e-6, 1e6)
#IsoO2._boundary_i.append([-0,IsoO2._conv_or_rad[-1]])
#IsoO2._boundary_i.insert(0, [1,IsoO2._conv_or_rad[0]])



#%% Structure Plots

# figsize, [lwr,lwc], ticksize, labelsize, legendfontsize, titlesize
PosterParams = [7.8, [2,4], 26, 24, 13, 20]

Plotter.Plot([IsoO1,MovO1,FreO1,FreO2], "p", "r", Params=PosterParams)
Plotter.Plot([IsoO1,MovO1,FreO1,FreO2], "k", "r", Params=PosterParams)
Plotter.Plot([IsoO1,MovO1,FreO1,FreO2], "T", "r", Params=PosterParams)

#Plotter.Plot([IsoO1,MovO1,FreO1], "p", "r", Params=PosterParams)
#Plotter.Plot([IsoO1,MovO1,FreO1], "k", "r", Params=PosterParams)
#Plotter.Plot([IsoO1,MovO1,FreO1], "T", "r", Params=PosterParams)

#Plotter.Plot([IsoO1,MovO1,FreO1,IsoO2], "p", "r", Params=PosterParams)
#Plotter.Plot([IsoO1,MovO1,FreO1,IsoO2], "k", "r", Params=PosterParams)
#Plotter.Plot([IsoO1,MovO1,FreO1,IsoO2], "T", "r", Params=PosterParams)



#%% Evolution Simulations

start_time = time.time()

# mc, me, rp, pd, kd, Td, L_guess_ini, rc, acc, t_tot, n_i, n_t
EIO1 = EIO.EvolverIsoO(12.07, 1.53, 40., 1e-10, 6e-5, 115, 1e-6, 3.6e9, 1e-10, 1.8e6, 1e4, 101)
EIO2 = EIO.EvolverIsoO(12.07, 1.53, 40., 1e-10, 1e-4, 115, 1e-6, 3.6e9, 1e-10, 3e6, 1e4, 101)
EIO3 = EIO.EvolverIsoO(12.07, 1.53, 40., 1e-10, 3e-4, 115, 1e-6, 3.6e9, 1e-10, 9e6, 1e4, 101)

end_time = time.time()
print("Runtime =", (end_time - start_time), "seconds")



#%% Evolution Plots

PosterParams = [7.8, 2, 26, 24, 13, 20]

Plotter.Plot([EIO1,EIO2,EIO3], "M", "t_pres", Params=PosterParams, Evo=True)
Plotter.Plot([EIO1,EIO2,EIO3], "L", "t_pres", Params=PosterParams, Evo=True)
Plotter.Plot([EIO2], "T_t", "r_t", Params=PosterParams, Evo=True)

#Plotter.Plot([EIO2], "phi_t", "r_t", Params=PosterParams, Evo=True)

















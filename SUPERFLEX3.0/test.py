# -*- coding: utf-8 -*-
"""
Created on Wed Aug 22 14:52:07 2018

@author: allen
"""

import numpy as np
import matplotlib.pyplot as plt

from superflex import Superflex
from Weigfun import Weigfun


#foring data
data=np.loadtxt('forcingWark.txt')
rainfall = data[:,4]
Qo = data[:,3]
evap = data[:,5]
#Q = Superflex(rainfall,evap, 0.1,0.1,0.1,0.2,0.2,0.2,300,0.5,0.4,0.2,0.3,100,0.2,0.1,
#          300,0.9,30,0.4,0.6,0.2,50,0.3,0.1,0.02,0.06,0.1,0.5,0.6,15,0.2,0.6,0.5,0.4,0.3,2.2,
#          1,1,1,1,1,1,1,1,1,1,1,1,1,1)
#option_i, option_u, option_f, option_s, option_c, res_i, res_r, res_u, res_f,res_s, res_c,
#             
        
#Calibration 
# Ce, K_Qq_FR, K_Qb_FR, m_E_FR, alpha_Qq_FR,
#              Beta_Qq_UR, Smax_UR, Beta_E_UR, SiniFr_UR, K_Qb_UR, mu_Qq_UR,
#              Smax_IR, m_QE_IR, K_Qq_RR, Smax_CR, Umax_uCR, Smin_uCR, Beta_Qq_uCR,
#              mu_Qq_uCR,K_Qb_uCR, Sevmax_CR, Beta_E_CR, Beta_Qq_sCR, K_Qb_sCR, K_Qd_sCR,
#              K_Qq_SR, m_E_SR, alpha_Qq_SR, P_ED_max,m_P_ED, D_S, D_I, D_F, D_R,
#              D_C, Tlag,
#                   Ce, K_Qq_FR,K_Qb_FR,m_E_FR ,alpha_Qq_FR, Beta_Qq_UR,Smax_UR
Par_min = np.array([0.01,0.01,  0.01,    0.01,      0.1,         0.1,      50,
                #   Beta_E_UR, SiniFr_UR, K_Qb_UR,mu_Qq_UR, Smax_IR, m_QE_IR
                    0.1,        0.01,      0.01,    0.1,     100,      0.1,
                #   K_Qq_RR, Smax_CR, Umax_uCR, Smin_uCR, Beta_Qq_uCR, mu_Qq_uCR
                    0.01,      200,     0.5,    10,           0.1,      0.1,
                #   K_Qb_uCR, Sevmax_CR, Beta_E_CR, Beta_Qq_sCR, K_Qb_sCR, K_Qd_sCR   
                    0.01,      10,         0.1,       0.1,        0.01,     0.01,
                #   K_Qq_SR, m_E_SR, alpha_Qq_SR, P_ED_max, m_P_ED, D_S, D_I, D_F
                    0.01,      0.1,   0.1,        10,       0.1,    0.1, 0.1, 0.1,
                #   D_R, D_C, T_lag
                    0.1, 0.1,   1])
#                   Ce,   K_Qq_FR, K_Qb_FR, m_E_FR ,alpha_Qq_FR, Beta_Qq_UR,Smax_UR
Par_max = np.array([0.5,   0.2,    0.2,    0.9,      0.9,         0.9,      500,
                #   Beta_E_UR, SiniFr_UR, K_Qb_UR,mu_Qq_UR, Smax_IR, m_QE_IR
                    0.9,        0.1,      0.1,    0.9,     500,      0.9,
                #   K_Qq_RR, Smax_CR, Umax_uCR, Smin_uCR, Beta_Qq_uCR, mu_Qq_uCR
                    0.1,      800,     0.9,     50,           0.9,      0.9,
                #   K_Qb_uCR, Sevmax_CR, Beta_E_CR, Beta_Qq_sCR, K_Qb_sCR, K_Qd_sCR   
                    0.1,       80,         0.9,       0.9,        0.1,     0.1,
                #   K_Qq_SR, m_E_SR, alpha_Qq_SR, P_ED_max, m_P_ED, D_S, D_I, D_F
                    0.1,      0.9,   0.9,        10000,       0.9,  0.9, 0.9, 0.9,
                #   D_R, D_C, T_lag
                    0.9, 0.9,   10])
n = 2000
obj_Par = []
for i in range(n):
    Par = Par_min + np.random.rand(len(Par_min))*(Par_max-Par_min)
    opt_ir = np.random.randint(1,3)
    opt_sr = np.random.randint(1,3)
    opt_fr = np.random.randint(1,3)
    opt_cr = np.random.randint(1,3)
    opt_ur = np.random.randint(1,5)
    lag_r = np.random.randint(0,2)
    lag_f = np.random.randint(0,2)
    lag_s = np.random.randint(0,2)
    res = np.random.randint(0,2,6)
    integer = np.append(np.array([lag_r,lag_f,lag_s,opt_ir,opt_ur,opt_fr,opt_sr,opt_cr]),res)
    Par = np.append(Par,integer)
    Obj,Qm = Superflex(rainfall,evap,Par,Qo)
    if Obj>.6 and Obj is not np.nan:
        obj_Par.append(Par)
        obj_Par[-1].append(Obj)
    #obj_Par.append(Obj)
    print("Processing %s/%s" %(i,n))
#Par = Par_min + np.random.rand(len(Par_min))*(Par_max-Par_min)
#opt_ir = np.random.randint(1,3)
#opt_sr = np.random.randint(1,3)
#opt_fr = np.random.randint(1,3)
#opt_cr = np.random.randint(1,3)
#opt_ur = np.random.randint(1,5)
#lag_r = 0#np.random.randint(0,1)
#lag_f = 0#np.random.randint(0,1)
#lag_s = 0#np.random.randint(0,1)
#res = np.random.randint(0,2,6)
#integer = np.append(np.array([lag_r,lag_f,lag_s,opt_ir,opt_ur,opt_fr,opt_sr,opt_cr]),res)
#Par = np.append(Par,integer)
#Par = np.append(Par,integer)
#obj,Qm = Superflex(rainfall,evap,Par,Qo)
#plt.figure()
#plt.plot(np.arange(0,len(Qo)),Qo,label='observed')
#plt.plot(np.arange(0,len(Qm)),Qm, label='modeled')
#plt.legend()
#plt.show()
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 25 11:13:37 2018

@author: allen
"""

import numpy as np
import matplotlib.pyplot as plt
from constfuns import *

def FUN_Mi(sEn_i,args):
    dt = args[0]
    P_i=args[1]
    Smax_i=args[2]
    m_i=args[3]
    option_i=args[4]
    Epot_i = args[5]
    if sEn_i<0:
        print('sEn_i<0')
#    global dt, P_i, Smax_i, m_i, option_i
    if option_i == 1:
        return dt*(P_i-P_i*power(sEn_i/Smax_i,m_i)-Epot_i*kinetics(sEn_i/Smax_i,m_i))
    elif option_i ==2:
        return dt*(P_i-P_i*rhfun(sEn_i/Smax_i,m_i)-Epot_i*kinetics(sEn_i/Smax_i,m_i))
    
    
def FUN_Mu(sEn_u, args):
    dt = args[0]
    P_u = args[1]
    Beta = args[2]
    mu = args[3]
    Gamma = args[4]
    Epot_u = args[5]
    Kb_u = args[6]
    option_u = args[7]
    Smax_u = args[8]
    
    if option_u ==1:
        return max(dt*(P_u-P_u*power(sEn_u/Smax_u, Beta)-Kb_u*(sEn_u/Smax_u)-Epot_u*kinetics(sEn_u/Smax_u,Gamma)),0)
    elif option_u ==2:
        return max(dt*(P_u-P_u*power(sEn_u/Smax_u, 1)-Kb_u*(sEn_u/Smax_u)-Epot_u*kinetics(sEn_u/Smax_u,Gamma)),0)
    elif option_u ==3:
        return max(dt*(P_u-P_u*kinetics(sEn_u/Smax_u,Beta)-Kb_u*(sEn_u/Smax_u)-Epot_u*kinetics(sEn_u/Smax_u,Gamma)),0)
    elif option_u ==4:
        return max(dt*(P_u-P_u*mlc(sEn_u/Smax_u,Beta,mu)-Kb_u*(sEn_u/Smax_u)-Epot_u*kinetics(sEn_u/Smax_u,Gamma)),0)

def FUN_MF(sEn_f,arg):
    dt = arg[0]
    P_f = arg[1]
    Kq_f = arg[2]
    Kb_f = arg[3]
    m_f = arg[4]
    alpha_f = arg[5]
    Epot_f = arg[6]
    D_FD = arg[7]
    option_f = arg[8]
    res_s = arg[9]
    res_u = arg[10]
    res_cr = arg[11]
    sSt_f = arg[12]
    
    if (res_u==0 and res_cr==0 and res_s==0):   #then UR, CR, SR are off
        if option_f ==1:
            return sEn_f-sSt_f-dt*(P_f-Kq_f*power(sEn_f,alpha_f)-Epot_f*Tfun(sEn_f,m_f))
        elif option_f ==2:
            return sEn_f-sSt_f-dt*(P_f-Kq_f*power(sEn_f,1)-Epot_f*Tfun(sEn_f,m_f))
    elif (res_u==0 and res_cr==1 and res_s==0):
        if option_f==1:
            return sEn_f-sSt_f-dt*(P_f-Kq_f*power(sEn_f,alpha_f))
        if option_f==2:
            return sEn_f-sSt_f-dt*(P_f-Kq_f*power(sEn_f,1))
    elif (res_u==1 and res_cr==0 and res_s==1):
        if option_f==1:
            return sEn_f-sSt_f-dt*(P_f-Kq_f*power(sEn_f,alpha_f)-Kb_f*sEn_f)
        if option_f==2:
            return sEn_f-sSt_f-dt*(P_f-Kq_f*power(sEn_f,1)-Kb_f*sEn_f)
    elif (res_u==1 and res_cr==1 and res_s==1):
        if option_f==1:
            return sEn_f-sSt_f-dt*(P_f-Kq_f*power(sEn_f,alpha_f)-Kb_f*sEn_f)
        if option_f==2:
            return sEn_f-sSt_f-dt*(P_f-Kq_f*power(sEn_f,1)-Kb_f*sEn_f)
    elif (res_u==1 and res_cr==0 and res_s==0):
        if option_f==1:
            return sEn_f-sSt_f-dt*(P_f-Kq_f*power(sEn_f,alpha_f))
        if option_f==2:
            return sEn_f-sSt_f-dt*(P_f-Kq_f*power(sEn_f,1))
    elif (res_u==1 and res_cr==1 and res_s==0):
        if option_f==1:
            return sEn_f-sSt_f-dt*(P_f-Kq_f*power(sEn_f,alpha_f))
        if option_f==2:
            return sEn_f-sSt_f-dt*(P_f-Kq_f*power(sEn_f,1))
    elif (res_u==0 and res_cr==0 and res_s==1):
        if option_f ==1:
            return sEn_f-sSt_f-dt*(P_f-Kq_f*power(sEn_f,alpha_f)-Kb_f*sEn_f-(1-D_FD)*Epot_f*Tfun(sEn_f,m_f))
        elif option_f ==2:
            return sEn_f-sSt_f-dt*(P_f-Kq_f*power(sEn_f,1)-Kb_f*sEn_f-(1-D_FD)*Epot_f*Tfun(sEn_f,m_f))
    elif (res_u==0 and res_cr==1 and res_s==1):
        if option_f==1:
            return sEn_f-sSt_f-dt*(P_f-Kq_f*power(sEn_f,alpha_f)-Kb_f*sEn_f)
        if option_f==2:
            return sEn_f-sSt_f-dt*(P_f-Kq_f*power(sEn_f,1)-Kb_f*sEn_f)  

    
    




def combined(sEn_ucr,arg):
    
    sSt_ucr = arg[0]
    sSt_scr = arg[1]
    P_ucr = arg[2]
    S_max_cr = arg[3]
    S_ev_max_cr = arg[4]
    U_max_ucr=arg[5]
    Epot_cr = arg[6]
    Beta_cr = arg[7]
    Beta_ucr=arg[8]
    mu_ucr=arg[9]
    Kb_ucr=arg[10]
    dt = arg[11]
    option_c = arg[11]
    S_min_ucr=arg[12]
    Beta_scr = arg[13]
    Kb_scr = arg[14]
    Kd_scr=arg[15]
    D_C = arg[16]
    
    sEn_ucr = max(sEn_ucr, S_min_ucr)
    sEn_scr = 0

    
    theta = sEn_ucr/(S_max_cr-sEn_scr)
    if option_c==1:
        Qq_ucr = P_ucr*power(theta/U_max_ucr,Beta_ucr)
    elif option_c==2:
        Qq_ucr = P_ucr*mlc(theta/U_max_ucr,Beta_ucr,mu_ucr)
    Qb_ucr = Kb_ucr*sEn_ucr
    if Qq_ucr+Qb_ucr > sEn_ucr-sSt_ucr:
        Qb_ucr= min(Qb_ucr,sEn_ucr - sSt_ucr) if Qb_ucr<sEn_ucr-sSt_ucr else max(0,sEn_ucr - sSt_ucr)
        Qq_ucr= max(0, sEn_ucr - sSt_ucr-Qb_ucr)

    P_scr = Qq_ucr*D_C+Qb_ucr
    sEn_scr = sSt_scr + P_scr
    Qq_scr = P_scr*power(sEn_scr/(S_max_cr-S_min_ucr),Beta_scr)
    Qb_scr = Kb_scr*sEn_scr
    Qd_scr = Kd_scr*sEn_scr
    if Qq_scr+Qb_scr+Qd_scr > sEn_scr-sSt_scr:
        Qb_scr = min(Qb_scr, sEn_scr-sSt_scr)
        Qd_scr = min(Qd_scr, sEn_scr-sSt_scr-Qb_scr) if sEn_scr-sSt_scr-Qb_scr>0 else 0
        Qq_scr = max(0,sEn_scr-sSt_scr-Qb_scr-Qd_scr) if sEn_scr-sSt_scr-Qb_scr-Qd_scr>0 else 0
    elif Qb_scr+Qd_scr>sEn_scr-sSt_scr:
        Qq_scr=0
        Qb_scr = min(Qb_scr, sEn_scr-sSt_scr)
        Qd_scr = max(0, sEn_scr-sSt_scr-Qb_scr)
    elif Qb_scr > sEn_scr-sSt_scr:
        Qq_scr=0
        Qd_scr=0
        Qb_scr = max(0, sEn_scr-sSt_scr)
#    print(Qb_ucr, Qq_ucr, Qq_scr, Qb_scr, Qd_scr)
    
    S_ev_cr = min(S_ev_max_cr, S_max_cr)
    S_ev_ucr = min(sEn_ucr,sEn_ucr+sEn_scr-S_ev_cr)
    S_ev_scr = max(0, sEn_scr-S_ev_cr)
    S_E_cr = S_ev_ucr+S_ev_scr
    S_E_cr_bar = S_E_cr/S_ev_cr
    
    E_cr = Epot_cr*kinetics(S_E_cr_bar,Beta_cr)
    
    E_ucr = min(E_cr*S_ev_scr/S_E_cr,sEn_ucr-sSt_ucr-Qb_ucr-Qq_ucr) if sEn_ucr-sSt_ucr-Qb_ucr-Qq_ucr>0 else 0
    E_scr = min(E_cr*S_ev_ucr/S_E_cr,sEn_scr-sSt_scr-Qb_scr-Qd_scr-Qq_scr) if sEn_scr-sSt_scr-Qb_scr-Qd_scr-Qq_scr>0 else 0
    Upright = dt*(P_scr-Qq_scr-Qb_scr-Qd_scr-E_scr)/(1-theta)
    Doright = dt*(P_ucr-Qb_ucr-Qq_ucr-E_ucr)-theta*(sEn_scr-sSt_scr)
    
    return sEn_ucr-sSt_ucr-Doright
    
    
def FUN_MI(sEn_i,arg):
    dt = arg[0]
    P_i=arg[1]
    Smax_i=arg[2]
    m_i=arg[3]
    option_i=arg[4]
    Epot_i = arg[5]
    sSt_i = arg[6]
    if sEn_i<0:
        print('sEn_i<0')
#    global dt, P_i, Smax_i, m_i, option_i
    if option_i == 1:
        return sEn_i-sSt_i-dt*(P_i-P_i*power(sEn_i/Smax_i,m_i)-Epot_i*kinetics(sEn_i/Smax_i,m_i))
    elif option_i ==2:
        return sEn_i-sSt_i-dt*(P_i-P_i*rhfun(sEn_i/Smax_i,m_i)-Epot_i*kinetics(sEn_i/Smax_i,m_i))
    
def FUN_MU(sEn_u,arg):
    dt = arg[0]
    P_u = arg[1]
    Beta = arg[2]
    mu = arg[3]
    Gamma = arg[4]
    Epot_u = arg[5]
    Kb_u = arg[6]
    option_u = arg[7]
    Smax_u = arg[8]
    sSt_u = arg[9]
    if option_u ==1:
        return sEn_u-sSt_u-dt*(P_u-P_u*power(sEn_u/Smax_u, Beta)-Kb_u*sEn_u-Epot_u*kinetics(sEn_u/Smax_u,Gamma))
    elif option_u ==2:
        return sEn_u-sSt_u-dt*(P_u-P_u*power(sEn_u/Smax_u, 1)-Kb_u*sEn_u-Epot_u*kinetics(sEn_u/Smax_u,Gamma))
    elif option_u ==3:
        return sEn_u-sSt_u-dt*(P_u-P_u*kinetics(sEn_u/Smax_u,Beta)-Kb_u*sEn_u-Epot_u*kinetics(sEn_u/Smax_u,Gamma))
    elif option_u ==4:
        return sEn_u-sSt_u-dt*(P_u-P_u*mlc(sEn_u/Smax_u,Beta,mu)-Kb_u*sEn_u-Epot_u*kinetics(sEn_u/Smax_u,Gamma))
    
def FUN_MS(sEn_s, arg):
    sSt_s = arg[0]
    dt = arg[1]
    P_s = arg[2]
    alpha_s = arg[3]
    Kq_s = arg[4]
    m_s = arg[5]
    D_F = arg[6]
    option_s = arg[7]
    E_pot_s = arg[8]
    res_c = arg[9]
    res_u = arg[10]
    res_f = arg[11]
#    print(option_s)
    if option_s ==1:
        Qq_SR = Kq_s*power(sEn_s, alpha_s)
    elif option_s ==2:
        Qq_SR = Kq_s*power(sEn_s, 1)
    if res_c==0 and res_u==0 and res_f==1:
        E_S = E_pot_s * D_F * Tfun(sEn_s, m_s)
    if res_c==1 or res_u==1:
        E_S = 0
    if res_c==0 and res_u==0 and res_f==0:
        E_S = E_pot_s *  Tfun(sEn_s, m_s)
    
    return sEn_s-sSt_s-P_s+Qq_SR+E_S

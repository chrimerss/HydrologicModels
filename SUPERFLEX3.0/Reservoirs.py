# -*- coding: utf-8 -*-
"""
Created on Sat Sep 15 10:49:05 2018

@author: allen
"""

import numpy as np
import matplotlib.pyplot as plt
from functions import *
from constfuns import power, kinetics, Tfun, rhfun, mlc
from scipy.optimize import fsolve

def InterceptionReservoir(dt, P_I, E_pot_I, option_I, S_I, Smax_I, m_I):
    sSt_I = S_I
    arg = [dt, P_I, Smax_I, m_I, option_I, E_pot_I, sSt_I]
    sEn_I = fsolve(FUN_MI, x0=sSt_I+P_I, args=(arg))
    E_I = E_pot_I*kinetics(sEn_I/Smax_I,m_I)
    
    if option_I==1:
        Q_IR = P_I*power(sEn_I/Smax_I, m_I)
    elif option_I == 2:
        Q_IR = P_I*rhfun(sEn_I/Smax_I, m_I)
    Q_IR = Q_IR+sEn_I-Smax_I if sEn_I>Smax_I else Q_IR
    sEn_I = min(Smax_I, sEn_I)
    if P_I-(sEn_I-sSt_I+Q_IR+E_I) >1e-3:
        print(P_I, sEn_I, sSt_i, Q_IR, E_I)
        raise ValueError("water is not balanced yet in IR")
    return sEn_I, Q_IR, E_I

def RiparianReservoir(dt, P_R, K_R, S_R):
    sSt_R = S_R
    sEn_R = (P_R*dt+sSt_R)/(1+K_R*dt)
    Q_RR = K_R *sEn_R
    
    if P_R-(sEn_R-sSt_R+Q_RR)>1e-3:
        print(P_R, sEn_R, sSt_R, Q_RR)
        raise ValueError("water is not balanced yet in RR")
    
    return sEn_R, Q_RR

def UnsaturatedReservoir(dt, P_U, E_pot_U, Beta, mu, Gamma, Kb_U, option_U, Smax_U, S_U ):
    sSt_U =S_U
    arg=[dt, P_U, Beta, mu, Gamma, E_pot_U, Kb_U, option_U, Smax_U, sSt_U]
    sEn_U = fsolve(FUN_MU, x0=sSt_U+P_U, args=(arg))
    if option_U ==1:
        Qq_UR = P_U*power(sEn_U/Smax_U, Beta)
    elif option_U ==2:
        Qq_UR = P_U*power(sEn_U/Smax_U, 1)
    elif option_U ==3:
        Qq_UR = P_U*kinetics(sEn_U/Smax_U,Beta)
    elif option_U ==4:
        Qq_UR = P_U * mlc(sEn_U/Smax_U,Beta,mu)
    E_U = E_pot_U*kinetics(sEn_U/Smax_U,Gamma)
    Qb_UR = Kb_U*sEn_U
#    print(P_U,sEn_U,sSt_U, Qq_UR)
    if P_U-(sEn_U-sSt_U+Qq_UR+Qb_UR+E_U)>1e-3:
        print(P_U, sEn_U, sSt_U, Qq_UR, Qb_UR, E_U)
        raise ValueError("water is not balanced yet in UR", P_U-(sEn_U-sSt_U+Qq_UR+Qb_UR))
    return sEn_U, Qq_UR, Qb_UR, E_U

def FastReservoir(dt, P_F, Kq_F, Kb_F, m_F, alpha_F, E_pot_F, D_F, option_F, res_u, res_c, S_F, res_s):
    sSt_F = S_F
    arg = [dt,P_F, Kq_F, Kb_F, m_F, alpha_F, E_pot_F, D_F, option_F,res_s, res_u, res_c, sSt_F ]
    sEn_F = fsolve(FUN_MF, x0=sSt_F+P_F, args=(arg))
    if sEn_F == 0:
        if res_c==0 and res_u==0 and res_s==0:
            E_F = E_pot_F * Tfun(sEn_F, m_F)
        elif res_u==0 and res_c==0 and res_s==1:
            E_F = E_pot_F* (1-D_F) * Tfun(sEn_F, m_F)
        elif res_u ==1 or res_c==1:
            E_F = 0
        Q_to_flow = sSt_F+P_F
        if option_F==1:
            Qq_FR = Kq_F*power(Q_to_flow,m_F)
        elif option_F==2:
            Qq_FR = Kq_F*power(Q_to_flow,1)
        if res_s==1:        
            Qb_FR = Q_to_flow-Qq_FR
        else:
            Qq_FR = Q_to_flow
            Qb_FR=0
        E_F = max(Q_to_flow-Qq_FR-Qb_FR,0)
    elif sEn_F != 0:
#            print(S_F[i])
        if res_c==0 and res_u==0 and res_s==0:
            E_F = E_pot_F * Tfun(sEn_F, m_F)
        elif res_u==0 and res_c==0 and res_s==1:
            E_F = E_pot_F* (1-D_F) * Tfun(sEn_F, m_F)
        elif res_u ==1 or res_c==1:
            E_F = 0
        if option_F == 1:
            Qq_FR = Kq_F*power(sEn_F, alpha_F)
        elif option_F ==2:
            Qq_FR = Kq_F*power(sEn_F, 1)
        if res_s==1:
            Qb_FR = Kb_F * sEn_F
        else:
            Qb_FR = 0
    if P_F-(sEn_F-sSt_F+Qb_FR+Qq_FR+E_F)>1e-3:
        print(P_F,sEn_F, sSt_F, Qb_FR, Qq_FR, E_F)
        raise ValueError("water is not balanced yet in FR")
        
    return sEn_F, Qb_FR, Qq_FR, E_F

def CombinedReservoir( P_C, S_UCR, S_SCR,E_pot_C, S_max_CR, S_evmax, U_max_UCR, Beta_CR, Beta_UCR, mu_UCR, Kb_UCR,
                      dt, option_C, S_min_UCR, Beta_SCR, Kb_SCR, Kd_SCR, D_C):
    sSt_UCR = max(S_UCR, S_min_UCR)
    sSt_SCR = S_SCR
    arg = [sSt_UCR,sSt_SCR,P_C,S_max_CR,S_evmax,U_max_UCR,E_pot_C, Beta_CR,
                   Beta_UCR,mu_UCR,Kb_UCR,dt,option_C,S_min_UCR,Beta_SCR,
                   Kb_SCR,Kd_SCR,D_C]
    sEn_UCR = fsolve(combined, x0=sSt_UCR+P_C, args=(arg))
    U_UCR = sEn_UCR/(S_max_CR- sEn_UCR)
    if option_C ==1:
        Qq_UCR = P_C*power(U_UCR/U_max_UCR,Beta_UCR)
    elif option_C ==2:
        Qq_UCR = P_C*mlc(U_UCR/U_max_UCR,Beta_UCR,mu_UCR)
    Qb_UCR = Kb_UCR*sEn_UCR
    flow_ucr = P_C-(sEn_UCR-sSt_UCR)
    if Qq_UCR+Qb_UCR > flow_ucr:
        Qb_UCR= min(Qb_UCR,flow_ucr) if Qb_UCR<flow_ucr else max(0,flow_ucr)
        Qq_UCR= max(0, flow_ucr-Qb_UCR)
    
    P_SCR = D_C* Qq_UCR + Qb_UCR
    sEn_SCR = P_SCR+sSt_SCR
    Qq_SCR = P_SCR*power(sEn_SCR/(S_max_CR-S_min_UCR), Beta_SCR)
    Qb_SCR = Kb_SCR*sEn_SCR
    Qd_SCR = Kd_SCR*sEn_SCR
    if Qq_SCR+Qb_SCR+Qd_SCR > sEn_SCR-sSt_SCR:
        Qb_SCR = min(Qb_SCR, sEn_SCR-sSt_SCR)
        Qd_SCR = min(Qd_SCR, sEn_SCR-sSt_SCR-Qb_SCR) if sEn_SCR-sSt_SCR-Qb_SCR>0 else 0
        Qq_SCR = max(0,sEn_SCR-sSt_SCR-Qb_SCR-Qd_SCR) if sEn_SCR-sSt_SCR-Qb_SCR-Qd_SCR>0 else 0
    elif Qb_SCR+Qd_SCR>sEn_SCR-sSt_SCR:
        Qq_SCR=0
        Qb_SCR = min(Qb_SCR, sEn_SCR-sSt_SCR)
        Qd_SCR = max(0, sEn_SCR-sSt_SCR-Qb_SCR)
    elif Qb_SCR > sEn_SCR-sSt_SCR:
        Qq_SCR=0
        Qd_SCR=0
        Qb_SCR = max(0, sEn_SCR-sSt_SCR)
    S_ev_CR = min(S_evmax, S_max_CR)
    S_ev_UCR = min(sEn_UCR,sEn_UCR+sEn_SCR-S_ev_CR) if sEn_UCR+sEn_SCR-S_ev_CR>0 else sEn_UCR
    S_ev_SCR = max(0, sEn_SCR-S_ev_CR)
    S_E_CR = S_ev_UCR+S_ev_SCR
    S_E_CR_bar = S_E_CR/S_ev_CR
    E_CR = E_pot_C*kinetics(S_E_CR_bar,Beta_CR)
    E_UCR = min(E_CR*S_ev_SCR/S_E_CR,sEn_UCR-sSt_UCR-Qb_UCR-Qq_UCR) if sEn_UCR-sSt_UCR-Qb_UCR-Qq_UCR>0 else 0
    E_SCR = min(E_CR*S_ev_UCR/S_E_CR,sEn_SCR-sSt_SCR-Qb_SCR-Qd_SCR-Qq_SCR) if sEn_SCR-sSt_SCR-Qb_SCR-Qd_SCR-Qq_SCR>0 else 0
    
    sEn_SCR -= Qq_SCR+Qd_SCR+Qb_SCR+E_SCR

    
    if P_C-(sEn_UCR-sSt_UCR+Qq_UCR+Qb_UCR+E_UCR)>1e-3:
        sEn_UCR += P_C-(sEn_UCR-sSt_UCR+Qq_UCR+Qb_UCR+E_UCR)
#    print(P_SCR-(sEn_SCR-sSt_SCR+Qq_SCR+Qb_SCR+Qd_SCR+E_SCR))
#    print(P_C-(sEn_UCR-sSt_UCR+Qq_UCR+Qb_UCR+E_UCR))
    sEn_C = sEn_UCR+sEn_SCR
    sSt_C = sSt_UCR+sSt_SCR
    E_C = E_UCR+E_SCR
    P_sl = (1-D_C)*Qq_UCR+Qq_SCR
    Q_CR = Qb_SCR+Qd_SCR
    
    if P_C-(sEn_C-sSt_C+P_sl+Q_CR+E_C)>1e-3:
        print(P_C, sEn_C, sSt_C, P_sl, Q_CR, E_C)
        raise ValueError("water is not balanced yet in CR")
        
    return sEn_UCR, sEn_SCR, E_C, Q_CR, P_sl

def SlowReservoir(dt, P_S, S_R, Kq_S, alpha_S, E_pot_S, D_F, m_S, res_c,option_S, res_u, res_f):
    sSt_S = S_R
    arg = [sSt_S, dt, P_S, alpha_S, Kq_S, m_S, D_F, option_S, E_pot_S,res_c, res_u, res_f]
    sEn_S = fsolve(FUN_MS, x0=sSt_S+P_S, args=(arg))
    if option_S ==1:
        Qq_SR = Kq_S*power(sEn_S, alpha_S)
    elif option_S ==2:
        Qq_SR = Kq_S*power(sEn_S, 1)
    if res_c==0 and res_u==0 and res_f==1:
        E_S = E_pot_S * D_F * Tfun(sEn_S, m_S)
    if res_c==1 or res_u==1:
        E_S = 0
    if res_c==0 and res_u==0 and res_f==0:
        E_S = E_pot_S * Tfun(sEn_S, m_S)
    if P_S-(Qq_SR+E_S+sEn_S-sSt_S)>1e-3:
        print(P_S, sEn_S, sSt_S, Qq_SR, E_S)
        raise ValueError("water is not balanced yet in SR")
    
    return sEn_S, Qq_SR, E_S

    
    
    

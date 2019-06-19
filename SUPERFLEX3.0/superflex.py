# -*- coding: utf-8 -*-
"""
Created on Sat Sep 15 11:54:15 2018

@author: allen
"""

import numpy as np
import matplotlib.pyplot as plt
from Reservoirs import InterceptionReservoir, RiparianReservoir, UnsaturatedReservoir, FastReservoir, SlowReservoir, CombinedReservoir
from Weigfun import Weigfun
from lagfun import tWeights
def Superflex(rainfall, evap, Par,Qo, dt=1):
    Ce = Par[0] 
    K_Qq_FR=Par[1]
    K_Qb_FR=Par[2]
    m_E_FR=Par[3]
    alpha_Qq_FR=Par[4]
    Beta_Qq_UR=Par[5]
    Smax_UR=Par[6]
    Beta_E_UR=Par[7]
    SiniFr_UR=Par[8]
    K_Qb_UR=Par[9]
    mu_Qq_UR=Par[10]
    Smax_IR=Par[11]
    m_QE_IR=Par[12]
    K_Qq_RR=Par[13]
    Smax_CR=Par[14] 
    Umax_uCR=Par[15]
    Smin_uCR=Par[16]
    Beta_Qq_uCR=Par[17]
    mu_Qq_uCR=Par[18]
    K_Qb_uCR=Par[19]
    Sevmax_CR=Par[20]
    Beta_E_CR=Par[21]
    Beta_Qq_sCR=Par[22]
    K_Qb_sCR=Par[23]
    K_Qd_sCR=Par[24]
    K_Qq_SR=Par[25]
    m_E_SR=Par[26]
    alpha_Qq_SR=Par[27]
    m_P_ED = Par[28]
    P_ED_max=Par[29]
    D_S=Par[30]
    D_I=Par[31]
    D_F=Par[32]
    D_R=Par[33]
    D_C=Par[34]
    Tlag=Par[35]
    Lag_RR=Par[36]
    Lag_FR=Par[37]
    Lag_SR=Par[38]
    option_i=Par[39] 
    option_u=Par[40]
    option_f=Par[41]
    option_s=Par[42]
    option_c=Par[43]
    res_i=Par[44]
    res_r=Par[45]
    res_u=Par[46] 
    res_f=Par[47]
    res_s=Par[48]
    res_c=Par[49] 
    
    P = rainfall
    evap = evap
    n=len(P)
    t = n
    #initialize

    # Interception reservoir
    Smax_I = Smax_IR
    m_i = m_QE_IR
    
    # Riparian reservoir
    K_r = K_Qq_RR
    lag_r = Lag_RR
    
    # Fast reservoir
    Kq_f = K_Qq_FR
    Kb_f = K_Qb_FR
    m_f = m_E_FR
    alpha_f = alpha_Qq_FR
    lag_f = Lag_FR
    
    #Slow reservoir
    Kq_s = K_Qq_SR
    m_s = m_E_SR
    alpha_s = alpha_Qq_SR
    lag_s = Lag_SR
    
    #unsaturated reservoir
    Beta = Beta_Qq_UR
    Smax_u = Smax_UR
    Gamma = Beta_E_UR
    Kb_u = K_Qb_UR
    mu = mu_Qq_UR
    
    # combined tank
       #UCR
    S_max_cr = Smax_CR
    U_max_ucr = Umax_uCR
    S_min_ucr = Smin_uCR
    Beta_ucr = Beta_Qq_uCR
    mu_ucr = mu_Qq_uCR
    Kb_ucr = K_Qb_uCR
       #SCR
    S_evmax = Sevmax_CR
    Beta_cr = Beta_E_CR
    Beta_scr = Beta_Qq_sCR
    Kb_scr = K_Qb_sCR
    Kd_scr = K_Qd_sCR
    #States
    S_I = np.zeros(n+1)
    S_U = np.zeros(n+1)
    S_U[0] = SiniFr_UR *Smax_u
    S_C = np.zeros(n+1)
    S_UCR = np.zeros(n+1)
    S_SCR = np.zeros(n+1)
    S_F = np.zeros(n+1)
    S_S = np.zeros(n+1)
    S_R = np.zeros(n+1)

    #Fluxes
    E_I = np.zeros(n+1)
    E_U = np.zeros(n+1)
    E_F = np.zeros(n+1)
    E_S = np.zeros(n+1)
    E_C = np.zeros(n+1)
    P_RL_vec = np.zeros(n+1)
    P_f_vec = np.zeros(n+1)
    P_sl_vec = np.zeros(n+1)
    Qtotal = np.zeros(n)
    for i in range(t):
#        i=1
        P_i = D_I*P[i]
        Q_ID = P[i]-P_i
        E_pot = evap[i]
        if res_i:
            S_I[i+1], Q_IR, E_I[i+1]=InterceptionReservoir(dt, P_i, E_pot, option_i, S_I[i], Smax_I, m_i)
        else:
            S_I[i+1]=E_I[i+1]=0
            Q_IR = P_i
        
        P_RD = Q_IR + Q_ID
        P_RL = D_R * P_RD
        P_RL_vec[i] = P_RL
        
        if lag_r:
            P_RL_lag = tWeights(Tlag, P_RL_vec, i+1)
        else:
            P_RL_lag = P_RL_vec
            
        P_r = P_RL_lag[i]
        if res_r:
            S_R[i+1], Q_RR= RiparianReservoir(dt, P_r, K_r, S_R[i])
        else:
            S_R[i+1]=0
            Q_RR=P_r
        # checking point 1
        
        
        P_SD = P_RD-P_RL
        Q_SD = D_S * P_SD
        P_ED = P_SD - Q_SD
        Q_ED = max(0, P_ED-P_ED_max)
        Q_ED = Q_ED/(1+np.exp((i+1)/m_P_ED))
        P_u = P_ED - Q_ED
        
        if res_u:
            E_pot_u = E_pot-Ce*E_I[i] if res_i==1 else E_pot
            S_U[i+1],Qq_UR, Qb_UR, E_U[i+1]=UnsaturatedReservoir(dt, P_u, E_pot_u, Beta, mu, Gamma, Kb_u, option_u, Smax_u, S_U[i])
        else:
            S_U[i+1]=S_U[i]
            E_U[i+1]=Qq_UR=0
            Qb_UR=P_u
        Q_UR = D_F*Qq_UR+Qb_UR
        P_f_vec[i] = Q_ED+Q_SD+(1-D_F)*Qq_UR

        if lag_f:
            P_f_lag = tWeights(Tlag,P_f_vec,i+1)
        else:
            P_f_lag = P_f_vec
        
        P_f = P_f_lag[i]
        if res_f:
            E_pot_f = E_pot - Ce*E_I[i]-Ce*E_U[i] 
            S_F[i+1], Qb_FR, Qq_FR, E_F[i+1]=FastReservoir(dt, P_f, Kq_f, Kb_f, m_f, alpha_f, E_pot_f, D_F, option_f, res_u, res_c, S_F[i], res_s)
        else:
            S_F[i+1]=Qb_FR=E_F[i+1]=0
            Qq_FR = P_f
        
        P_c = Q_UR
        if res_c:
            S_UCR[i] = max(S_UCR[i], S_min_ucr)
            S_C[i] = max(S_UCR[i],S_C[i])
            E_pot_c = E_pot-Ce*E_I[i]-Ce*E_U[i]
            S_UCR[i+1],S_SCR[i+1], E_C[i+1], Q_CR, P_sl = CombinedReservoir(P_c, S_UCR[i], S_SCR[i],E_pot_c, S_max_cr, S_evmax, U_max_ucr, Beta_cr, Beta_ucr, mu_ucr, Kb_ucr,
                      dt, option_c, S_min_ucr, Beta_scr, Kb_scr, Kd_scr, D_C)
            S_C[i+1] = S_UCR[i+1]+S_SCR[i+1]
        else:
            S_C[i+1]=E_C[i+1]=Q_CR=0
            P_sl = P_c
        
        P_sl_vec[i] = P_sl
        
        if lag_s:
            P_sl_lag = tWeights(Tlag,P_sl_vec,i+1)
        else:
            P_sl_lag = P_sl_vec
            
        P_s = P_sl_lag[i]+Qb_FR
        if res_s:
            E_pot_s = E_pot -Ce*E_F[i]- Ce*E_I[i]
            S_S[i+1],Qq_SR, E_S[i+1]= SlowReservoir(dt, P_s, S_S[i], Kq_s, alpha_s, E_pot_s, D_F, m_s, res_c,option_s, res_u, res_f)
        else:
            S_S[i+1]=E_S[i+1]=0
            Qq_SR = P_s
        
        Qtotal[i] = Qq_SR+Qq_FR+Q_CR+Q_RR
        storage_change = S_I[i+1]+S_R[i+1]+S_U[i+1]+S_C[i+1]+S_F[i+1]+S_S[i+1]-S_I[i]-S_R[i]-S_U[i]-S_C[i]-S_F[i]-S_S[i]
#        if P[i]-(storage_change+Qtotal[i]+E_I[i+1]+E_U[i+1]+E_C[i+1]+E_F[i+1]+E_S[i+1]) >1e-3:
#            raise ValueError("water doesnot balance in time step")
#        print("water balance at each step:",P[i]-(storage_change+Qtotal[i]+E_I[i+1]+E_U[i+1]+E_C[i+1]+E_F[i+1]+E_S[i+1]))
#    Weigths=Weigfun(Tlag)
#    Qm = np.convolve(Qtotal,Weigths)
    Qm=Qtotal
#    if len(Qm) != len(Qo):
#        print( len(Qtotal),len(Qm),len(Qo))
#        raise ValueError("The length of modeled and observed doesnot match")
#    Qm = Qtotal
    # NS objective function
    ind=np.where(Qo>=0)
    QoAv=np.mean(Qo[ind])
    ErrUp=sum((Qm[ind]-Qo[ind])**2)
    ErrDo=sum((Qo[ind]-QoAv)**2)
    Obj=1-ErrUp/ErrDo
    print(Obj)
    Initial_U = SiniFr_UR *Smax_u
    Initial_C = S_min_ucr if res_c==1 else 0
    Initial_S =  Initial_U+Initial_C
#    print("Initial Storage is",Initial_S)
#    print("-----------------------------------------------------------------")
#    print("The total discharge is %.2f"%(Qtotal.sum()))
#    print("The storage value of SI is %.2f"%(S_I[-1]))
#    print("The storage value of SR is %.2f"%(S_R[-1]))
#    print("The storage value of SC is %.2f"%(S_C[-1]))
#    print("The storage value of SU is %.2f"%(S_U[-1]))
#    print("The storage value of SF is %.2f"%(S_F[-1]))
#    print("The storage value of SS is %.2f"%(S_S[-1]))
    Total_E = E_I.sum()+E_U.sum()+E_C.sum()+E_F.sum()+E_S.sum()
#    print("Total Evaporation is %.2f"%(Total_E))
#    print("potential evaporation is ",evap.sum())
#    print("E_C is %.2f"%(E_C.sum()))
    Total_S = S_I[-1]+S_R[-1]+S_C[-1]+S_U[-1]+S_F[-1]+S_S[-1]
#    print("water balance %.2f"%(Total_S+Total_E+Qtotal.sum()-P.sum()-Initial_S))
#    print(res_i, res_r, res_u, res_f, res_s, res_c)
    return Obj,Qm

    
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
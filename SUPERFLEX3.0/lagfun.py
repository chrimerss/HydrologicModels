# -*- coding: utf-8 -*-
"""
Created on Wed Aug 22 11:18:10 2018

@author: allen
"""

import numpy as np
import matplotlib.pyplot as plt

def tWeights(Nlag, Qin, n):
    wc = int(np.ceil(Nlag))
    wf = int(np.floor(Nlag))
    LagArr = np.zeros(wc)
    Qout = np.zeros(n)
    
    m =  2.0/Nlag**2
    
    w = np.asarray([(2*i+1)*m for i in range(wc)])
    if wf<wc:
        w[wc-1] = (m*Nlag**2-m*wf**2)/2.0
    
    w = w/sum(w)
    
    for i in range(n):
        LagArr = LagArr +Qin[i] *w
        Qout[i] = LagArr[0]
        
        LagArr[:wc-1] = LagArr[1:wc]
        LagArr[wc-1] = 0
    if np.any(Qout==np.nan):
        print("Yes")
    return Qout


# -*- coding: utf-8 -*-
"""
Created on Wed Aug 22 14:39:17 2018

@author: allen
"""

import numpy as np
import matplotlib.pyplot as plt

def power(x,m):
#    if x<0:
#        x=0
    return x**m

def kinetics(x,m):
    return x*(1+m)/(x+m)

def Tfun(x,m):
    return 1-np.exp(-x/m)

def rhfun(x,m):
    return 1-(1-x)*(1+m)/(1-x+m)

def mlc(x,m,lab):
    return (1+np.exp(-m*(1-lab)))*(np.exp(-m*x)-1)/ (1+np.exp(-m*(x-lab)))/(np.exp(-m)-1)  
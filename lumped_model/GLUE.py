# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 10:52:50 2018

@author: allen
"""

import numpy as np
import matplotlib.pyplot as plt
from HBVMod import HBVMod
import pandas as pd

forcing=np.genfromtxt('Forcing.txt',  dtype=float, autostrip=True)[:,3:6]

# Monte Carlo Simulation
n=10000
ParMinn = np.array([0,   0.2,  40,    .5,   .001,   0,     .01,  .0001])
ParMaxn = np.array([8,    1,  800,   4,    .3,     10,    .1,   .01])
Sin= np.array([0,  100,  0,  5  ])

A=np.zeros((n,9))
n_feasible = 0
Par=np.zeros(len(ParMinn))
Q = []
for i in range(1,n): 
    Rnum=np.random.rand(len(ParMinn)) #generate a vector of random number
    Par = Rnum * (ParMaxn - ParMinn) + ParMinn
# calculate the random parameter set
    Obj,Qm = HBVMod(Par, forcing, Sin, 'false') #call the model
    print('proceeding {} %'.format(i/n*100,2))
    if Obj>.6:
        A[n_feasible,0:8]= Par
        A[n_feasible,8]=Obj
        n_feasible = n_feasible + 1
        Q.append(Qm)
A = A[:n_feasible]
data_list=list(A)
# GLUE to pick out feasible parameter
#step one: define likelihood function 
#step two: set threshold to generate behavioal realization
#step three: rescale likelihood
#data = pd.DataFrame(columns=['Index','Q','L','P'])
## cumulative sum
#Q_max=[]
#Q_min=[]
#for t in range(len(forcing[:,0])):
#    Q_vec = []
#    data = pd.DataFrame(columns=['Index','Q','L','P'])
#    for j in range(len(data_list)):
#        Qm = Q[j][t]
#        Q_vec.append(Qm)
#    Q_vec_sort = np.sort(Q_vec)
#    index = np.argsort(Q_vec)
#    Likelihood = A[:,8][index]
#    P = np.cumsum(Likelihood)/sum(Likelihood)
#    data.Index = index
#    data.Q = Q_vec_sort
#    data.L = Likelihood
#    data.P = P
#    index_first = data.Index[data['P']<0.05]
#    index_last = data.Index[data.P>0.95]
#    data= data[data.P >0.05]
#    data = data[data.P<0.95]
#    Q_max.append(max(data.Q))
#    Q_min.append(min(data.Q))
#    print(t)
#length = len(forcing[:,0])
#plt.figure()
#plt.plot(np.arange(length),Q_max, color='r')
#plt.plot(np.arange(length),Q_min,color='r')
#plt.plot(np.arange(length),forcing[:,1])
#plt.show()

class GLUE:
    def __init__(self,Qm,Qo, calibrated):
        '''
        Qm: the modeled discharge data; list(t, Qt)
        Qo: the observed discharge data
        calibrated: nd-array like; including parameters and likelihood (the last)
        '''
        self.Qm = Qm
        self.Qo = Qo
        self.Qmin = []
        self.Qmax = []
        self.A = calibrated
        self._score=0
        
    def uncertainty(self,percentile, method='normal'):
        if method=='normal':
            for t in range(len(self.Qo)):
                Q_vec = []
                data = pd.DataFrame(columns=['Index','Q','L','P'])
                for j in range(len(self.A[:,-1])):
                    _Qm = self.Qm[j][t]
                    Q_vec.append(_Qm)
                Q_vec_sort = np.sort(Q_vec)
                index = np.argsort(Q_vec)
                Likelihood = self.A[:,-1][index]
                P = np.cumsum(Likelihood)/sum(Likelihood)
                data.Index = index
                data.Q = Q_vec_sort
                data.L = Likelihood
                data.P = P
  #              index_head = data.Index[data.P< percentile]
 #               index_tail = data.Index[data.P> percentile]
                data= data[data.P >percentile]
                data = data[data.P<1-percentile]
                self.Qmax.append(max(data.Q))
                self.Qmin.append(min(data.Q))
            
    
    def score(self, method='CR'):
 #       if not isinstance(self.Qo, np.array):
 #           raise ValueError("Qo should be in right format")
        if not isinstance(self.Qmin,list):
            raise ValueError('Qmin is not defined')
        if not isinstance(self.Qmax,list):
            raise ValueError('Qmax is not defined')
        if len(self.Qmin)==0 or len(self.Qmax)==0:
            raise ValueError('Qmin or Qmax havenot defined yet')
        if method == 'CR':
            #for t in range(len(self.Qo)):
            score = [1 if np.all(self.Qo[t]<self.Qmax[t] and self.Qo[t]> self.Qmin[t]) else 0 for t in range(len(self.Qo))]
            score = np.array(score)
            self._score = (score/len(self.Qo)).sum()
            return self._score
            

    def plot(self):
        plt.figure()
        plt.plot(np.arange(len(self.Qo)), self.Qmax, color='r', label='maximum')
        plt.plot(np.arange(len(self.Qo)), self.Qmin, color='r', label='minmum')
        plt.plot(np.arange(len(self.Qo)),self.Qo, label='observed')
        plt.legend()
        plt.text(1,1000,'score:%s'%(self._score))
        plt.show()

glue=GLUE(Q, forcing[:,1],A)
glue.uncertainty(0.05)
glue.score()
glue.plot()
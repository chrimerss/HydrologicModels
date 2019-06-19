import numpy as np
import matplotlib as mpl
from HBVMod import HBVMod

forcing=np.genfromtxt('Forcing.txt',  dtype=float, autostrip=True)

          #      Imax Ce Sumax beta Pmax   Tlag   Kf  Ks
ParMinn = np.array([0,   0.2,  40,    .5,   .001,   0,     .01,  .0001])
ParMaxn = np.array([8,    1,  800,   4,    .3,     10,    .1,   .01])
Sin= np.array([0,  100,  0,  5  ])

forcing= forcing[:,3:6]


# GLUE
nmax=5000
A=np.zeros((nmax,9))
Q=[]
n_feasible = 0
Par=np.zeros(len(ParMinn))
for n in range(1,nmax): 
    for i in range(len(ParMinn)):
        Rnum=np.random.rand(1) #generate a vector of random number
        Par[i] = Rnum * (ParMaxn[i] - ParMinn[i]) + ParMinn[i]
# calculate the random parameter set
    Obj,Qm = HBVMod(Par, forcing, Sin, 'false') #call the model
    print('proceeding {} %'.format(n/nmax*100,2))
    if Obj>.6:
        A[n_feasible,0:8]= Par
        A[n_feasible,8]=Obj
        n_feasible = n_feasible + 1
        Q.append(Qm)
    

np.savetxt('MC2.txt',A[0:n_feasible,:], delimiter =',')


#find the optimum
#find the optimal parameter set

Obj,Qm =HBVMod(Par,forcing,Sin, hydrograph='TRUE')

# GLUE
t = np.arange(0,1725,5)
perf = A[:len(Q),8]
perf = np.sort(perf)
weight_perf = [ 1/len(perf) * i for i in range(len(perf))]
dis = Q
# rescale weights
perf = (perf[:] - min(perf) )/ (max(perf)-min(perf))
for time in t:
    discharge = dis[:][time]
    P = [discharge , perf]
    discharge_rank = np.sort(discharge)
    weight_discharge = np.array([ 1/len(discharge_rank) * i for i in range(len(discharge_rank))])
    discharge = discharge_rank[(weight_discharge>0.05) & (weight_discharge<0.95)]
    dis=dis[dis[:][time]==discharge]
    
Q.index(dis)











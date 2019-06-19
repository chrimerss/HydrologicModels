import numpy       as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from HBVMod import HBVMod


forcing=np.genfromtxt('wark_data/forcingWark.txt',  dtype=float, autostrip=True)
# we consider it as in the hillslope so Pmax change to D
                  #      Imax Ce Sumax beta  D Tlag   Kf  Ks 
#Par = np.array([6.15, 0.68, 89.15, 1.85, 0.09, 1.10, 0.1, 0.008])
#Par_min = np.array([3,0.3,100,0.9 ,0.2,1 ,0.93, 0.01])
#Par_max = np.array([9,0.8,400,0.99,0.6,5 ,0.99, 0.04])
#Consider as normal case
Par = np.array([6.15, 0.68, 89.15, 1.85, 0.09, 1.10, 0.1, 0.008])
Par_min = np.array([3,0.3,100,0.9 ,0.01,1 ,0.93, 0.01])
Par_max = np.array([9,0.8,400,0.99,0.3,5 ,0.99, 0.04])
Par_difference = Par_max-Par_min
num = 10000
objective = []
Par_objective = []
for i in range(num):
    rand = np.random.rand(8)
    Par = Par_difference*rand +Par_min
    Qm = HBVMod(Par, forcing[:,3:6])
    Qo = forcing[:,3]
    ind=np.where(Qo>=0)
    QoAv=np.mean(Qo[ind])
    ErrUp=sum((Qm[ind]-Qo[ind])**2)
    ErrDo=sum((Qo[ind]-QoAv)**2)
    Obj=1-ErrUp/ErrDo
    if Obj>0.6:
        objective.append(Obj)
        Par_objective.append(Par)
    print("Proceeding %s/10000"%(i))
np.corrcoef([Qm,Qo])

ind = objective.index(max(objective))
Par = Par_objective[ind]
print(max(objective))
print(Par)
Qm = HBVMod(Par, forcing[:,3:6])
plt.plot(range(0,len(Qo)),Qo)
plt.plot(range(0,len(Qm)),Qm)
plt.legend(['Observed','Modeled'])
plt.xlabel('Time Series')
plt.ylabel('Discharge mm/day')
plt.title('Comparison of Measured discharge and modeling')
plt.show()
import numpy       as np
import matplotlib.pyplot as plt
from HBVMod import HBVMod
plt.style.use('ggplot')
forcing=np.genfromtxt('Forcing.txt',  dtype=float, autostrip=True)

Sin= np.array([0,  100,  0,  5  ])

forcing= forcing[:,3:6]

A=np.genfromtxt('MC2.txt',  dtype=float, autostrip=True, delimiter=',')


#find the optimum
#find the optimal parameter set
OptPar = A[A[:,8]==max(A[:,8])]

plt.figure(1)
Obj,Qm =HBVMod(OptPar[0][0:8],forcing,Sin, hydrograph='TRUE')

plt.figure(2)
plt.subplot(421)
plt.plot(A[:,0],A[:,8],'o')
plt.xlabel('Interception')
plt.ylabel('Objective')

plt.subplot(422)
plt.plot(A[:,1],A[:,8],'o')
plt.xlabel('Ce')
plt.ylabel('Objective')

plt.subplot(423)
plt.plot(A[:,2],A[:,8],'o')
plt.xlabel('Sumax')
plt.ylabel('Obejective')


plt.subplot(424)
plt.plot(A[:,3],A[:,8],'o')
plt.xlabel('beta')
plt.ylabel('Objective')


plt.subplot(425)
plt.plot(A[:,4],A[:,8],'o')
plt.xlabel('Pmax')
plt.ylabel('Objective')

plt.subplot(426)
plt.plot(A[:,5],A[:,8],'o')
plt.xlabel('Tlag')
plt.ylabel('Objective')

plt.subplot(427)
plt.plot(A[:,6],A[:,8],'o')
plt.xlabel('Kf')
plt.ylabel('Obejective')

plt.subplot(428)
plt.plot(A[:,7],A[:,8],'o')
plt.xlabel('Ks')
plt.ylabel('Obejective')

plt.show()

Obj,Qm=HBVMod(OptPar[0][0:8],forcing,Sin, hydrograph='s')
#  Box Plot

plt.figure
plt.subplot(421)
plt.boxplot(A[:,0])
plt.title('Interception')

plt.subplot(422)
plt.boxplot(A[:,1])
plt.title('Ce')

plt.subplot(423)
plt.boxplot(A[:,2])
plt.title('Sumax')


plt.subplot(424)
plt.boxplot(A[:,3])
plt.title('beta')


plt.subplot(425)
plt.boxplot(A[:,4])
plt.title('Pmax')


plt.subplot(426)
plt.boxplot(A[:,5])
plt.title('Tlag')

plt.subplot(427)
plt.boxplot(A[:,6])
plt.title('Kf')

plt.subplot(428)
plt.boxplot(A[:,7])
plt.title('Ks')

plt.show()
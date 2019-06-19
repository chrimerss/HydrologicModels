import numpy       as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from FLEXtopo import FLEXtopo


forcing=np.genfromtxt('wark_data/forcingWark.txt',  dtype=float, autostrip=True)

                  #      Imax Ce Sumax beta Pmax    Kf  
ParPlateau   = np.array([3.2, 0.50, 17.40, 0.95, 1.76, 0.91])   
                  #      Imax Ce Sumax beta D     Kf  
ParHillslope = np.array([3.25, 0.50, 321.99, 0.99, 0.4,0.97])
                  #      Imax Ce Sumax beta Cmax     Kf  
ParWetland   = np.array([9.94, 0.50, 53.25, 0.70, 0.65, 0.45])
              # Ks Tlag
ParCatchment = np.array([0.0281, 2.21])

#landscape percentages
landscape_per=np.array([hillslope_per,wetland_per,plateau_per])
num = 1000
ParPlateau_min = np.array([2.8,0.3, 10,0.85,1,0.85])
ParPlateau_max = np.array([9,  0.8, 50,0.98,2,0.93])
ParPlateau_difference = ParPlateau_max-ParPlateau_min
ParHillslope_min = np.array([3,0.3,100,0.9 ,0.2,0.93])
ParHillslope_max = np.array([9,0.8,400,0.99,0.6,0.99])
ParHillslope_difference = ParHillslope_max-ParHillslope_min
ParWetland_min = np.array([8 ,0.3,25,0.5,0.4,0.3])
ParWetland_max = np.array([12,0.8,75,0.9,0.8,0.6])
ParCatchment_min = np.array([0.01,1])
ParCatchment_max = np.array([0.04,5])
ParCatchment_difference = ParCatchment_max - ParCatchment_min
ParWetland_difference = ParWetland_max-ParWetland_min
Par_difference = [ParPlateau_difference, ParHillslope_difference, ParWetland_difference]
Par_min = [ParPlateau_min, ParHillslope_min, ParWetland_min]
Par = [ParPlateau , ParHillslope , ParWetland]
Par_optimal=[]
objective_Par = []
objective_Cat = []
for i in range(num):
    rand_catchment = np.random.rand(2)
    ParCatchment = rand_catchment*ParCatchment_difference + ParCatchment_min
    for j in range(3):
        rand_par = np.random.rand(6)
        Par[j] = Par_difference[j] *rand_par + Par_min[j]
    Qm = FLEXtopo(Par[0], Par[1], Par[2], ParCatchment, forcing[:,3:6], landscape_per )
    Qo = forcing[:,3]
    ind=np.where(Qo>=0)
    QoAv=np.mean(Qo[ind])
    ErrUp=sum((Qm[ind]-Qo[ind])**2)
    ErrDo=sum((Qo[ind]-QoAv)**2)
    Obj=1-ErrUp/ErrDo
    if Obj > 0.6:
        Par_optimal.append(Par)
        objective_Par.append(Obj)
        objective_Cat.append(ParCatchment)
    print("proceeding %s /%s"%(i,num))

    # Calculate objective
ind = objective_Par.index(max(objective_Par))
print(max(objective_Par))
Par_obj = Par_optimal[ind]
Cat_obj = objective_Cat[ind]
Qm = FLEXtopo(Par_obj[0], Par_obj[1], Par_obj[2], Cat_obj, forcing[:,3:6], landscape_per )
plt.plot(range(0,len(Qo)),Qo)
plt.plot(range(0,len(Qm)),Qm)
plt.legend(['Observed','Modeled'])
plt.xlabel('Time Series')
plt.ylabel('Discharge mm/day')
plt.title('Comparison of Measured discharge and modeling')
plt.show()
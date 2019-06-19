import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


DEM=np.genfromtxt('E:/Learning Materials/TUD/hydrologic modeling/week 6/AS 6/python_files/exercise_files/wark_data/dem.asc', 
                  dtype=float, autostrip=True)
slope=np.genfromtxt('E:/Learning Materials/TUD/hydrologic modeling/week 6/AS 6/python_files/exercise_files/wark_data/slope.asc', 
                  dtype=float, autostrip=True)
hand=np.genfromtxt('E:/Learning Materials/TUD/hydrologic modeling/week 6/AS 6/python_files/exercise_files/wark_data/HAND.asc', 
                  dtype=float, autostrip=True)
basin=np.genfromtxt('E:/Learning Materials/TUD/hydrologic modeling/week 6/AS 6/python_files/exercise_files/wark_data/basin.asc', 
                  dtype=float, autostrip=True)

#plot DEM
plt.figure(1)
DEM[DEM==-9999]=np.nan
plt.imshow(DEM, cmap='hsv')
plt.colorbar()
 
#plot HAND
plt.figure(2)
hand[hand==-9999]=np.nan
plt.imshow(hand,cmap='hsv')
plt.colorbar()


#make landscape classification
hillslope = (np.array(slope)) > 11
plateau = (np.array(hand) > 5) & (np.array(slope)  <11)
wetland = (np.array(hand) )<= 5 &( np.array(hand)>0)
basin = (np.array(basin))>0

hillslope_per = float(np.sum(hillslope))/float(np.sum(basin))
wetland_per = float(np.sum(wetland))/float(np.sum(basin))
plateau_per = float(np.sum(plateau))/float(np.sum(basin))

landscapes=np.zeros((118,220))
landscapes[plateau]=1
landscapes[hillslope]=2
landscapes[wetland]=3

#plot landscapes
cmap = mpl.colors.ListedColormap(['white', 'red', 'green', 'blue'])
bounds=[0,1,2,3,4]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

plt.figure(3)
plt.imshow(landscapes, cmap=cmap, norm=norm)
plt.colorbar()
plt.show()






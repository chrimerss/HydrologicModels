import numpy       as np
from Weigfun import Weigfun
from plateau import plateau
from hillslope import hillslope
from wetland import wetland

def FLEXtopo( ParPlateau, ParHillslope, ParWetland, ParCatchment,forcing, landscapes):


	#parameters and constants
	Tlag=ParCatchment[1]
	Ks=ParCatchment[0]
	dt=1
	tmax=len(forcing[:,0])

	#initialize states
	States_plateau=np.zeros((tmax,3))
	States_hillslope=np.zeros((tmax,3))
	States_wetland=np.zeros((tmax,3))
	Ss=np.zeros((tmax,1))

	#initialize fluxes
	Fluxes_plateau=np.zeros((tmax,4))
	Fluxes_hillslope=np.zeros((tmax,4))
	Fluxes_wetland=np.zeros((tmax,3))
	Qsdt=np.zeros(tmax)
	Qtotdt=np.zeros(tmax)

	#
	#loop over time
	for t in range(0,tmax):
		
		#plateau
		Fluxes_plateau, States_plateau=plateau( t, ParPlateau, forcing, Fluxes_plateau, States_plateau )
		#hillslope
		Fluxes_hillslope, States_hillslope=hillslope( t, ParHillslope, forcing, Fluxes_hillslope, States_hillslope )

		#wetland
		Fluxes_wetland, States_wetland, Ss=wetland( t, ParWetland, forcing, Fluxes_wetland, States_wetland, Ss, landscapes[2] )

		# Slow Reservoir
		Ss[t]=Ss[t]+ Fluxes_plateau[t,3] * landscapes[0] + Fluxes_hillslope[t,3] * landscapes[1] 
		Qsdt= dt*Ks*Ss[t] 
		Ss[t]=Ss[t]-min(Qsdt,Ss[t])
		if t<tmax-1:
			Ss[t+1]=Ss[t]

		Qtotdt[t]=Qsdt + Fluxes_plateau[t,2] * landscapes[0] + Fluxes_hillslope[t,2] * landscapes[1] +  Fluxes_wetland[t,2] * landscapes[2] 




	# Offset Q

	Weigths=Weigfun(Tlag)
	
	Qm = np.convolve(Qtotdt,Weigths)
	Qm=Qm[0:tmax]
	
	return(Qm)



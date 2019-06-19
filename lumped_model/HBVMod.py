import numpy       as np
import matplotlib.pyplot as plt
from Weigfun import Weigfun

def HBVMod( Par,forcing,Sin, hydrograph):
	#HBVpareto Calculates values of 3 objective functions for HBV model

    Imax=Par[0]
    Ce=Par[1]
    Sumax=Par[2]
    beta=Par[3]
    Pmax= Par[4]
    Tlag= Par[5]
    Kf= Par[6]
    Ks= Par[7]

    
    Prec=forcing[:,0]
    Qo=forcing[:,1]
    Etp=forcing[:,2]


    tmax=len(Prec)
	
	# allocate Si, Su, Sf, Ss, Eidt, Eadt, Qtotdt
    Si=np.zeros(tmax)
    Su=np.zeros(tmax)
    Sf=np.zeros(tmax)
    Ss= np.zeros(tmax)
    Qtotdt=np.zeros(tmax)
    Eadt=np.zeros(tmax)
    Eidt=np.zeros(tmax)

	# initialize Si, Su, Sf, Ss
    Si[0]=Sin[0]
    Su[0]= 100
    Ss[0]= 5

    dt=1

	#
	# Model 1 SOF1
    for i in range(0,tmax):
        Pdt=Prec[i]*dt
        Epdt=Etp[i]*dt
	    # Interception Reservoir
        if Pdt>0:
            Si[i]=Si[i]+Pdt
            if Si[i]> Imax:
                Pedt = Si[i] - Imax
                Si[i] = Imax
            else:
                Si[i]=Si[i]
                Pedt = 0
            Eidt[i]= 0
        else:
		# Evaporation only when there is no rainfall
            Pedt= 0
            Eidt[i] = min(Si[i],Epdt)
            Si[i]= Si[i] - Eidt[i]
	   
        if i<tmax-1:
            Si[i+1]=Si[i]
	    
	    
	    # Unsaturated Reservoir
        if Pedt>0:
            rho=(Su[i]/Sumax)**beta
            Qiudt=(1-rho)*Pedt*dt # flux from Ir to Ur
            Su[i]=Su[i]+Qiudt
            Qufdt=rho*Pedt*dt  #flux from Su to Sf
        else:
            Qufdt=0
            Qiudt=0
	    
	    # Transpiration
        Epdt=max(0,Epdt-Eidt[i])
        Eadt[i]= Su[i]/Sumax/Ce *dt* Epdt
        Eadt[i]=min(Eadt[i], Su[i])  #
        Su[i]= Su[i] - Eadt[i]
	    # Percolation
        Qusdt=Pmax*(Su[i]/Sumax)*dt # Flux from Su to Ss
        Su[i]=Su[i] - Qusdt
        if i<tmax-1:
            Su[i+1]=Su[i]
	    
	    # Fast Reservoir
        Sf[i]=Qufdt + Sf[i]
        Qfdt= dt*Kf*Sf[i]
        Sf[i]= Sf[i] - Qfdt
        if i<tmax-1:
            Sf[i+1]=Sf[i]
	    
	    # Slow Reservoir
        Ss[i]= Ss[i] + Qusdt
        Qsdt= Ks* Ss[i] *dt
        Ss[i]= Ss[i] - Qsdt
        if i<tmax-1:
            Ss[i+1]=Ss[i]
	    
        Qtotdt[i]=Qsdt+Qfdt


	# Check Water Balance
    Sf=Si[-1]+Ss[-1]+Sf[-1]+Su[-1] #final storage
    Sin=sum(Sin) #initial storage
    WB=sum(Prec)-sum(Eidt)-sum(Eadt) -sum(Qtotdt)- Sf +Sin
	# Offset Q

    Weigths=Weigfun(Tlag)
	
    Qm = np.convolve(Qtotdt,Weigths)
    Qm=Qm[0:tmax]
    # Calculate objective
    ind=np.where(Qo>=0)
    QoAv=np.mean(Qo[ind])
    ErrUp=sum((Qm[ind]-Qo[ind])**2)
    ErrDo=sum((Qo[ind]-QoAv)**2)
    Obj=1-ErrUp/ErrDo
    if hydrograph == 'TRUE':
	## Plot
	# hour=1:tmax
        plt.figure
        plt.plot(range(0,len(Qo)),Qo)
        plt.plot(range(0,len(Qm)),Qm)
        plt.legend(['Observed','Modeled'])
        plt.xlabel("Time series")
        plt.ylabel("Discharge mm/day")
        plt.title("Comparison of measured discharge and modelling")
        plt.show()
#    print(np.corrcoef([Qo,Qm]))
    return(Obj,Qm)

	
	# leg['Qobs','Qmod']

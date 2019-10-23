import numpy as np
import math as math
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

#-------> FUNCTIONS <-------

def FSLE(timeseries,do= 1,num=70,r=1.1):
        
    GridGoals = np.zeros([num])
    GridGoals[0] = do
    FSLE = np.zeros([3,num])         #I would like to keep the second row as index of how many pairs compute that mean, you see what I mean?
    
    for i in range (0,num-1) :
        GridGoals[i+1] = GridGoals[i]*r
        
    Errors=[]
    for i in range(0,num) :
        n = 0
        Time = 0 
        for t in range(0,len(timeseries)) : 
            A = timeseries[t]
            for k in range(1,len(A[0][:])) : 
                if A[0][k]> GridGoals[i] and A[0][k-1]< GridGoals[i]: 
                    Time = Time + (np.abs(A[1][k]-A[1][0]))
                    n = n+1 
                    Errors.append(np.abs(A[1][k]-A[1][0]))
                    break
        if n != 0 :
            FSLE[0,i] = math.pow(Time/n,-1)*math.log(r)
            FSLE[1,i] = n
            Var = 0 
            for p in range(0,len(Errors)) :
                Var = Var + ((Time/n - Errors[p])*(Time/n - Errors[p]))/n    
                FSLE[2,i] = math.pow(Var,0.5) * math.pow(Time/n,-2)*math.log(r)
        else :
            print("Warning  :",i)
        while len(Errors) > 0 : Errors.pop()
    return GridGoals,FSLE

def Kd(timeseries,Dbeg=-1,Dend=3,nbin=50):
    MaxLen = 1
    for i in range(len(timeseries)):
        if timeseries[i].shape[1] > MaxLen:
            MaxLen  =  timeseries[i].shape[1]
    K = np.zeros((len(timeseries),MaxLen,2))
    
    for i in range(len(timeseries)):
        for j in range(len(timeseries[i][0])-1):
            K[i,j,0] = ((timeseries[i][0][j+1]*1000)**2-(timeseries[i][0][j]*1000)**2)/3600/2
            K[i,j,1] = (timeseries[i][0][j]+timeseries[i][0][j+1])/2
            
    Dbin = np.logspace(Dbeg,Dend,nbin)
    Kx = Dbin+Dend/nbin/2
    Kstat = np.zeros((3,nbin))
    for i in range(nbin-1):
        n = 0
        Ks = []
        for j in range(len(K)):
            for k in range(len(K[j])):
                if K[j,k,1]<Dbin[i+1] and K[j,k,1]>Dbin[i]:
                    Ks += [K[j,k,0]]
                    n += 1
        if n!=0:
            Kstat[0,i] = np.mean(Ks)
            Kstat[1,i] = n
            Kstat[2,i] = np.std(Ks)
                       
    return Kx,Kstat

#-------> MAIN <-------

DirecInputForward = "Data/ForwardDistances/"
DirecInputBackward = "Data/BackwardsDistances/"

PairAsy = []
PairSyFor =  []
PairSyBac =  []
ChancePairFor = []  
ChancePairBac = []

for i in range (0,33):
   A = np.genfromtxt(DirecInputBackward+'BDPair'+str(i)+".csv", delimiter=',')  #It read the file and copy in A
   Ab = np.genfromtxt(DirecInputForward+'FDPair'+str(i)+".csv", delimiter=',') 
   if len(A.shape) == 1 :
      PairAsy.append(Ab)
   else : 
      if A[0][-1] < 5 :
         PairSyFor.append(Ab) 
         PairSyBac.append(A)
      else : 
         ChancePairFor.append(Ab)
         ChancePairBac.append(A)

############################---FSLE---###############################

GridGoals, FSLE_Asy = FSLE(PairAsy)
GridGoals, FSLE_SyFor = FSLE(PairSyFor)
GridGoals, FSLE_SyBac = FSLE(PairSyBac)
GridGoals, FSLE_ChanceFor = FSLE(ChancePairFor)
GridGoals, FSLE_ChanceBac = FSLE(ChancePairBac)


xError = np.zeros(len(GridGoals)) #Error space (equal to zero)

#############################---Diffusivity---################################
Kx,KPA = Kd(PairAsy)
Kx,KPF = Kd(PairSyFor)
Kx,KPB = Kd(PairSyBac)
Kx,KCF = Kd(ChancePairFor)
Kx,KCB = Kd(ChancePairBac)

#-------> PLOT <-------
fig, ax1 = plt.subplots()

ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_xlabel('Distance (km)')
ax1.set_ylabel('FSLE')

ax1.errorbar(GridGoals, FSLE_SyFor[0,:],
            xerr=xError,
            yerr=FSLE_SyFor[2,:],
            fmt='-o', color=(0.5, 0.1, 0.4, 0.9))
ax1.errorbar(GridGoals, FSLE_SyBac[0,:],
            xerr=xError,
            yerr=FSLE_SyBac[2,:],
            fmt='-o', color=(0.5, 0.1, 0.4, 0.3))
ax1.errorbar(GridGoals, FSLE_Asy[0,:],
            xerr=xError,
            yerr=FSLE_Asy[2,:],
            fmt='-o', color=(0.1, 0.5, 0.2, 0.9))
ax1.errorbar(GridGoals, FSLE_ChanceFor[0,:],
            xerr=xError,
            yerr=FSLE_ChanceFor[2,:],
            fmt='-o', color=(0.9, 0.3, 0., 0.6))
ax1.errorbar(GridGoals, FSLE_ChanceBac[0,:],
            xerr=xError,
            yerr=FSLE_ChanceBac[2,:],
            fmt='-o', color=(0.9, 0.3, 0., 0.2))

ax2 = ax1.twinx() 


ax2.set_ylabel('Components')  
ax2.plot(GridGoals, FSLE_SyFor[1,:], color=(0.5, 0.1, 0.4, 0.9))
ax2.plot(GridGoals, FSLE_SyBac[1,:], color=(0.5, 0.1, 0.4, 0.3))
ax2.plot(GridGoals, FSLE_Asy[1,:], color=(0.1, 0.5, 0.2, 0.9))
ax2.plot(GridGoals, FSLE_ChanceFor[1,:], color=(0.9, 0.3, 0., 0.6))
ax2.plot(GridGoals, FSLE_ChanceBac[1,:], color=(0.9, 0.3, 0., 0.2))


fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()

#%%
 
fig, ax4 = plt.subplots()
ax4.set_yscale('log')
ax4.set_xscale('log')
ax4.set_ylabel('Diffusivity K [$m^2$/s]')
ax4.set_xlabel('Distance (km)')
ax4.plot(Kx,KPA[0],'o',color='brown',label='Asy')
ax4.fill_between(Kx,KPA[0]-KPA[2],KPA[0]+KPA[2],color='brown',alpha=0.3)
ax4.loglog(Kx,KPB[0],'<',color='gray',label='SyBack')
ax4.fill_between(Kx,KPB[0]-KPB[2],KPB[0]+KPB[2],color='gray',alpha=0.3)
ax4.loglog(Kx,KPF[0],'>',color='pink',label='SyForw')
ax4.fill_between(Kx,KPF[0]-KPF[2],KPF[0]+KPF[2],color='pink',alpha=0.3)
ax4.loglog(Kx,KCB[0],'<',color='blue',label='ChanceBack')
ax4.fill_between(Kx,KCB[0]-KCB[2],KCB[0]+KCB[2],color='blue',alpha=0.3)
ax4.loglog(Kx,KCF[0],'>',color='red',label='ChanceForw')
ax4.fill_between(Kx,KCF[0]-KCF[2],KCF[0]+KCF[2],color='red',alpha=0.3)

plt.legend()
plt.show()
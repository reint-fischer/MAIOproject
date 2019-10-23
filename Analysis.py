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

#############################---Diffusivity---################################
def Kd(timeseries):
    MaxLen = 1
    for i in range(len(timeseries)):
        if timeseries[i].shape[1] > MaxLen:
            MaxLen  =  timeseries[i].shape[1]
    K = np.zeros((len(timeseries),MaxLen,2))
    
    for i in range(len(timeseries)):
        for j in range(len(timeseries[i][0])-1):
            K[i,j,0] = (timeseries[i][0][j+1]**2-timeseries[i][0][j]**2)/3600/2
            K[i,j,1] = (timeseries[i][0][j]+timeseries[i][0][j+1])/2
    return K

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

#-------> PLOT <-------

fig, ax1 = plt.subplots()

ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_xlabel('Distance (km)')
ax1.set_ylabel('FSLE')

ax1.errorbar(GridGoals, FSLE_SyFor[0,:],
            xerr=xError,
            yerr=FSLE_SyFor[2,:],
            fmt='-o', color='pink')
ax1.errorbar(GridGoals, FSLE_SyBac[0,:],
            xerr=xError,
            yerr=FSLE_SyBac[2,:],
            fmt='-o', color='gray')
ax1.errorbar(GridGoals, FSLE_Asy[0,:],
            xerr=xError,
            yerr=FSLE_Asy[2,:],
            fmt='-o', color='brown')
ax1.errorbar(GridGoals, FSLE_ChanceFor[0,:],
            xerr=xError,
            yerr=FSLE_ChanceFor[2,:],
            fmt='-o', color='red')
ax1.errorbar(GridGoals, FSLE_ChanceBac[0,:],
            xerr=xError,
            yerr=FSLE_ChanceBac[2,:],
            fmt='-o', color='blue')

ax2 = ax1.twinx() 


ax2.set_ylabel('Components')  
ax2.plot(GridGoals, FSLE_SyFor[1,:], color='pink')
ax2.plot(GridGoals, FSLE_SyBac[1,:], color='gray')
ax2.plot(GridGoals, FSLE_Asy[1,:], color='brown')
ax2.plot(GridGoals, FSLE_ChanceFor[1,:], color='red')
ax2.plot(GridGoals, FSLE_ChanceBac[1,:], color='blue')


fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()


####Backward###
#f2 = plt.figure(2)
#ax2 = plt.axes()
#plt.title('Backward FSLE')
#color = 'tab:red'
#ax2.set_xlabel('Distance (km)')
##ax1.set_yscale('log')
#ax2.set_ylabel('FSLE', color=color)
#ax2.loglog(GridGoals/1000, FSLE[0,:,1], color=color)
#ax2.errorbar(GridGoals/1000, FSLE[0,:,1],
#            xerr=xError,
#            yerr=FSLE[2,:,1],
#            fmt='-o', color='red')
#ax2.tick_params(axis='y', labelcolor=color)
#
#ax22 = ax2.twinx()  # instantiate a second axes that shares the same x-axis
#
#color = 'tab:blue'
#ax22.set_ylabel('Components', color=color)  # we already handled the x-label with ax1
#ax22.plot(GridGoals/1000 , FSLE[1,:,1], color=color)
#ax22.tick_params(axis='y', labelcolor=color)
#
#f2.tight_layout()  # otherwise the right y-label is slightly clipped
#
####Comparison###
#f3 = plt.figure(3)
#ax3 = plt.axes()
#plt.title('Symmetry')
#color = 'tab:red'
#ax3.set_xlabel('Distance (km)')
##ax1.set_yscale('log')
#ax3.set_ylabel('FSLE')
#ax3.loglog(GridGoals/1000, FSLE[0,:,0], color=color)
#ax3.errorbar(GridGoals/1000, FSLE[0,:,0],
#            xerr=xError,
#            yerr=FSLE[2,:,0],
#            fmt='-o', color='red')
#
#color = 'tab:green'
#ax3.loglog(GridGoals/1000, FSLE[0,:,1], color=color)
#ax3.errorbar(GridGoals/1000, FSLE[0,:,1],
#            xerr=xError,
#            yerr=FSLE[2,:,1],
#            fmt='-o', color='green')
#
#

K = Kd(Pair) 
f4 = plt.figure(4,figsize=(10,8))
ax4 = plt.axes()
ax4.set_ylabel('Diffusivity K [$m^2$/s]')
ax4.set_xlabel('Distance (km)')
for i in K:  
    ax4.loglog(i[:,1],i[:,0],'o')
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
                    Time = Time + (np.abs(A[1][k]-A[1][0])/24)
                    n = n+1 
                    Errors.append(np.abs(A[1][k]-A[1][0])/24)
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
ax1.set_ylabel('FSLE [day^-1]')

ax1.plot(GridGoals, FSLE_SyFor[0,:], '-o', color='pink', label='SyForw')
ax1.fill_between(GridGoals, FSLE_SyFor[0,:]+(FSLE_SyFor[2,:]/2), FSLE_SyFor[0,:]-(FSLE_SyFor[2,:]/2), facecolor='pink', alpha=0.5)     

ax1.plot(GridGoals, FSLE_SyBac[0,:], '-o', color='gray', label='SyBack')
ax1.fill_between(GridGoals, FSLE_SyBac[0,:]+(FSLE_SyBac[2,:]/2), FSLE_SyBac[0,:]-(FSLE_SyBac[2,:]/2), facecolor='gray', alpha=0.5)   

ax1.plot(GridGoals, FSLE_Asy[0,:], '-o', color='brown', label='Asy')
ax1.fill_between(GridGoals, FSLE_Asy[0,:]+(FSLE_Asy[2,:]/2), FSLE_Asy[0,:]-(FSLE_Asy[2,:]/2), facecolor='brown', alpha=0.5)  

ax1.plot(GridGoals, FSLE_ChanceFor[0,:], '-o', color='red', label='ChanceForw')
ax1.fill_between(GridGoals, FSLE_ChanceFor[0,:]+(FSLE_ChanceFor[2,:]/2), FSLE_ChanceFor[0,:]-(FSLE_ChanceFor[2,:]/2), facecolor='red', alpha=0.5)  

ax1.plot(GridGoals, FSLE_ChanceBac[0,:], '-o', color='blue', label='ChanceBack')
ax1.fill_between(GridGoals, FSLE_ChanceBac[0,:]+(FSLE_ChanceBac[2,:]/2), FSLE_ChanceBac[0,:]-(FSLE_ChanceBac[2,:]/2), facecolor='blue', alpha=0.5)  

ax2 = ax1.twinx() 
ax2.set_ylabel('Components')  

ax2.plot(GridGoals, FSLE_SyFor[1,:], color='pink')
ax2.plot(GridGoals, FSLE_SyBac[1,:], color='gray')
ax2.plot(GridGoals, FSLE_Asy[1,:], color='brown')
ax2.plot(GridGoals, FSLE_ChanceFor[1,:], color='red')
ax2.plot(GridGoals, FSLE_ChanceBac[1,:], color='blue')

fig.tight_layout()  # otherwise the right y-label is slightly clipped
legend = ax1.legend(loc='upper right', shadow=True, fontsize='x-large')

plt.show()

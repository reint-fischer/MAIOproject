import numpy as np
import math as math
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

DirecInputForward = "Data/ForwardDistances/"
DirecInputBackward = "Data/BackwardsDistances/"
pairs = np.genfromtxt('Data/UnPair.txt', delimiter=',')



ChancePair = []       #-> Two list, one with the pairs and one with the chancepairs
Pair = []
MaxLenPair = 0        #-> For localizing the dataset which is the longest in the group
MaxLenChancePair = 0

ChancePairb = []       #-> Two list, one with the pairs and one with the chancepairs
Pairb = []
MaxLenPairb = 0        #-> For localizing the dataset which is the longest in the group
MaxLenChancePairb = 0

for i in range (0,33):
    Ab = np.genfromtxt(DirecInputBackward+'BDPair'+str(i)+".csv", delimiter=',')  #It read the file and copy in A
    A = np.genfromtxt(DirecInputForward+'FDPair'+str(i)+".csv", delimiter=',')  #It read the file and copy in A
    if len(Ab.shape) > 1:
        ChancePairb.append(Ab)
        ChancePair.append(A)
        Bb = len(Ab[0][:])
        B = len(A[0][:])
        if Bb > MaxLenChancePairb :
            MaxLenChancePairb = Bb
        if B > MaxLenChancePair :
            MaxLenChancePair = B
    else:
        Pairb.append(Ab)
        Bb = 1
        if Bb > MaxLenPairb :
            MaxLenPairb = Bb
        Pair.append(A)
        B = len(A[0][:])
        if B > MaxLenPair :
            MaxLenPair = B
        
    

############################---FSLE---################################

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

GridGoals,Fd = FSLE(ChancePair)
xError = np.zeros(len(GridGoals)) #Error space (equal to zero)


###Forward###
f1 = plt.figure(1)
ax1 = plt.axes()
plt.title('Forward FSLE')
color = 'tab:red'
ax1.set_xlabel('Distance (km)')
#ax1.set_yscale('log')
ax1.set_ylabel('FSLE', color=color)
ax1.loglog(FSLE(ChancePair)[0], FSLE(ChancePair)[1][0], color=color)
ax1.errorbar(FSLE(ChancePair)[0], FSLE(ChancePair)[1][0],
            xerr=xError,
            yerr=FSLE(ChancePair)[1][2],
            fmt='-o', color='red')
ax1.tick_params(axis='y', labelcolor=color)

ax12 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax12.set_ylabel('Components', color=color)  # we already handled the x-label with ax1
ax12.plot(FSLE(ChancePair)[0], FSLE(ChancePair)[1][1], color=color)
ax12.tick_params(axis='y', labelcolor=color)

f1.tight_layout()  # otherwise the right y-label is slightly clipped

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
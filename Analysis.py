import numpy as np
import math as math
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

DirecInputForward = "Data/ForwardDistances/"
DirecInputBackward = "Data/BackwardsDistances/"
pairs = np.genfromtxt('Data/UnPair.txt', delimiter=',')

def func(x, a, b, c):
   return a * math.pow(x, -2/b) + c   #THIS FOR THE FIT -> IT DOESN'T WORK


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


do = 1000 #[meters]
num = 70 #Number of resolution FSLE
r = 1.1 #Grid increment

GridGoals = np.zeros([num])
GridGoals[0] = do
FSLE = np.zeros([3,num,2])         #I would like to keep the second row as index of how many pairs compute that mean, you see what I mean?

for i in range (0,num-1) :
   GridGoals[i+1] = GridGoals[i]*r

Errors=[]
for i in range(0,num) :
   n = 0
   Time = 0 
   for t in range(0,len(ChancePair)) : 
      A = ChancePair[t]
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

#popt, pcov = curve_fit(func, GridGoals, FSLE[0,:])
#Model = np.zeros([num])
#for i in range (0,num):
#   Model[i] = popt[0] * math.pow(x, -2/popt[1]) + popt[2]
#print(popt)
   
Errorsb=[]
for i in range(0,num) :
   nb = 0
   Timeb = 0 
   for t in range(0,len(ChancePairb)) : 
      Ab = ChancePairb[t]
      for k in range(1,len(Ab[0][:])) : 
         if Ab[0][k]> GridGoals[i] and Ab[0][k-1]< GridGoals[i]: 
            Timeb = Timeb + (np.abs(Ab[1][k]-Ab[1][0]))
            nb = nb+1 
            Errorsb.append(np.abs(Ab[1][k]-Ab[1][0]))
            break
   if nb != 0 :
      FSLE[0,i,1] = math.pow(Timeb/nb,-1)*math.log(r)
      FSLE[1,i,1] = nb
      Varb = 0 
      for p in range(0,len(Errorsb)) :
         Varb = Varb + ((Timeb/nb - Errorsb[p])*(Timeb/nb - Errorsb[p]))/nb    
      FSLE[2,i,1] = math.pow(Varb,0.5) * math.pow(Timeb/nb,-2)*math.log(r)
   else :
      print("Warning  :",i)
   while len(Errorsb) > 0 : Errorsb.pop()

xError = np.zeros(num) #Error space (equal to zero)


###Forward###
f1 = plt.figure(1)
ax1 = plt.axes()
plt.title('Forward FSLE')
color = 'tab:red'
ax1.set_xlabel('Distance (km)')
#ax1.set_yscale('log')
ax1.set_ylabel('FSLE', color=color)
ax1.loglog(GridGoals/1000, FSLE[0,:,0], color=color)
ax1.errorbar(GridGoals/1000, FSLE[0,:,0],
            xerr=xError,
            yerr=FSLE[2,:,0],
            fmt='-o', color='red')
ax1.tick_params(axis='y', labelcolor=color)

ax12 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax12.set_ylabel('Components', color=color)  # we already handled the x-label with ax1
ax12.plot(GridGoals/1000 , FSLE[1,:,0], color=color)
ax12.tick_params(axis='y', labelcolor=color)

f1.tight_layout()  # otherwise the right y-label is slightly clipped

###Backward###
f2 = plt.figure(2)
ax2 = plt.axes()
plt.title('Backward FSLE')
color = 'tab:red'
ax2.set_xlabel('Distance (km)')
#ax1.set_yscale('log')
ax2.set_ylabel('FSLE', color=color)
ax2.loglog(GridGoals/1000, FSLE[0,:,1], color=color)
ax2.errorbar(GridGoals/1000, FSLE[0,:,1],
            xerr=xError,
            yerr=FSLE[2,:,1],
            fmt='-o', color='red')
ax2.tick_params(axis='y', labelcolor=color)

ax22 = ax2.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax22.set_ylabel('Components', color=color)  # we already handled the x-label with ax1
ax22.plot(GridGoals/1000 , FSLE[1,:,1], color=color)
ax22.tick_params(axis='y', labelcolor=color)

f2.tight_layout()  # otherwise the right y-label is slightly clipped

###Comparison###
f3 = plt.figure(3)
ax3 = plt.axes()
plt.title('Symmetry')
color = 'tab:red'
ax3.set_xlabel('Distance (km)')
#ax1.set_yscale('log')
ax3.set_ylabel('FSLE')
ax3.loglog(GridGoals/1000, FSLE[0,:,0], color=color)
ax3.errorbar(GridGoals/1000, FSLE[0,:,0],
            xerr=xError,
            yerr=FSLE[2,:,0],
            fmt='-o', color='red')

color = 'tab:green'
ax3.loglog(GridGoals/1000, FSLE[0,:,1], color=color)
ax3.errorbar(GridGoals/1000, FSLE[0,:,1],
            xerr=xError,
            yerr=FSLE[2,:,1],
            fmt='-o', color='green')


############################---Diffusivity---################################
K = np.zeros((len(Pair),MaxLenPair,2))
f4 = plt.figure(4,figsize=(10,8))
ax4 = plt.axes()
ax4.set_ylabel('Diffusivity K [$m^2$/s]')
ax4.set_xlabel('Distance (km)')
for i in range(len(Pair)):
    for j in range(len(Pair[i][0])-1):
        K[i,j,0] = (Pair[i][0][j+1]**2-Pair[i][0][j]**2)/3600/2
        K[i,j,1] = (Pair[i][0][j]+Pair[i][0][j+1])/2
    ax4.loglog(K[i,:,1]/1000,K[i,:,0],'o',label = 'floats {0} and {1}'.format(int(pairs[i,0]),int(pairs[i,1])))

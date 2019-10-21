import numpy as np
import math as math
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def func(x, a, b, c):
   return a * math.pow(x, -2/b) + c   #THIS FOR THE FIT -> IT DOESN'T WORK

DirecInputForward = "/home/giovanni/MAIO/Python/ForwardDistances/"
DirecInputBackward = "/home/giovanni/MAIO/Python/BackwardsDistances/"

ChancePair = []       #-> Two list, one with the pairs and one with the chancepairs
Pair = []
MaxLenPair = 0        #-> For localizing the dataset which is the longest in the group
MaxLenChancePair = 0

for i in range (0,33):
   A = np.genfromtxt(DirecInputForward+'FDPair'+str(i)+".csv", delimiter=',')  #It read the file and copy in A
   if min(A[0][:]) < 1 :                                                       #If the minimum value is below 1km (squared) it's a pair
      Pair.append(A)
      B = len(A[0][:])
      if B > MaxLenPair :
         MaxLenPair = B 
   if min(A[0][:]) > 1 :                                                       #If the minimum distance is above 1km (squared) it's a chancepair
      ChancePair.append(A)
      B = len(A[0][:])
      if B > MaxLenChancePair :
         MaxLenChancePair = B 

############################---FSLE---################################

do = 0.5     #[kilometers]
num = 50     #Number of resolution FSLE
r = 1.1      #Grid increment
GridGoals = np.zeros([num])
GridGoals[0] = do
FSLE = np.zeros([3,num])         #I would like to keep the second row as index of how many pairs compute that mean, you see what I mean?

for i in range (0,num-1) :
   GridGoals[i+1] = GridGoals[i]*r

Errors=[]
for i in range(0,num) :
   n = 0
   Time = 0 
   for t in range(0,len(ChancePair)) : 
      A = ChancePair[t]
      for k in range(1,len(A[0][:])) : 
         if A[0][k] > GridGoals[i]*GridGoals[i] and A[0][k-1] < GridGoals[i]*GridGoals[i]: 
            Time = Time + (A[1][k]-A[1][0])
            n = n+1 
            Errors.append(A[1][k]-A[1][0])
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

xError = np.zeros(num) #Error space (equal to zero)
fig, ax1 = plt.subplots()
color = 'tab:red'
ax1.set_xlabel('Distance (km)')
ax1.set_yscale('log')

ax1.set_ylabel('FSLE', color=color)
ax1.errorbar(GridGoals, FSLE[0,:],
            xerr=xError,
            yerr=FSLE[2,:],
            fmt='-o', color='red')
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('Components', color=color)  # we already handled the x-label with ax1
ax2.plot(GridGoals, FSLE[1,:], color=color)
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()


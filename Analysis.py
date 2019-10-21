import numpy as np
import math as math
import matplotlib.pyplot as plt

DirecInputForward = "Data/ForwardDistances/"
DirecInputBackward = "Data/BackwardsDistances/"

ChancePair = []
Pair = []
MaxLenPair = 0
MaxLenChancePair = 0

for i in range (0,32):
   A = np.genfromtxt(DirecInputForward+'FDPair'+str(i)+".csv", delimiter=',')/1000
   if min(A[0][:]) < 1 : 
      Pair.append(A)
      B = len(A[0][:])
      if B > MaxLenPair :
         MaxLenPair = B 
   if min(A[0][:]) > 1 : 
      ChancePair.append(A)
      B = len(A[0][:])
      if B > MaxLenChancePair :
         MaxLenChancePair = B 

############################---FSLE---################################

do = 1.5 #[meters]
num = 65 #Number of resolution FSLE
r = 1.1 #Grid increment
GridGoals = np.zeros([num])
GridGoals[0] = do
FSLE = np.zeros([2,num])  #I would like to keep the second row as index of how many pairs compute that mean, you see what I mean?
for i in range (0,num-1) :
   GridGoals[i+1] = GridGoals[i]*r

for i in range(0,num) :
   n = 0
   Time = 0 
   for t in range(0,len(ChancePair)) : 
      A = ChancePair[t]
      for k in range(1,len(A[0][:])) : 
         if A[0][k]> GridGoals[i] and A[0][k-1]< GridGoals[i]: 
            Time = Time + (A[1][k]-A[1][0])
            n = n+1 
            break
   if n != 0 :
      FSLE[0,i] = math.pow(Time/n,-1)*math.log(r)
      FSLE[1,i] = n
   else :
      print("Warning  :",i)


fig, ax1 = plt.subplots()
color = 'tab:red'
ax1.set_xlabel('Distance (km)')
#ax1.set_yscale('log')
ax1.set_ylabel('FSLE', color=color)
ax1.loglog(GridGoals, FSLE[0,:], color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('Components', color=color)  # we already handled the x-label with ax1
ax2.plot(GridGoals, FSLE[1,:], color=color)
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()

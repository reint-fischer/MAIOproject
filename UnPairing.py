import numpy as np

Pair = np.genfromtxt('Data/Pair.txt', delimiter=',')


k = 0 
UnPair = np.zeros([50,2])
for i in range (0,70):
    A = Pair[i]
    B = [A[1] ,A[0]]
    if not (B in Pair): print('Fuck it') #There's always another
    if not( A in UnPair) : 
       if not ( B in UnPair) : 
          UnPair[k] = A
          k= k+1

np.savetxt('Data/UnPair.txt',UnPair,delimiter=',')

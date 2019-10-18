import numpy as np

def SelectPair(D):
    Pair= np.zeros([100,2])
    k=0
    for i in range ( int(len(D[2,:])/4) ): 
        print( 100*(i*1.0)/(1.0*int(len(D[2,:])/4)), '%' )
        for t in range ( int(len(D[2,:])/4) ):
            if D[0,i*4]!= D[0,t*4] : 
                space_distance = (D[2,i*4]-D[2,t*4])*(D[2,i*4]-D[2,t*4]) + (D[1,i*4]-D[1,t*4])*(D[1,i*4]-D[1,t*4])
                time_distance = (D[3,i*4] - D[3,t*4])*(D[3,i*4] - D[3,t*4]) 
                if space_distance < 100 and time_distance < 36 :
                   if not( [D[0,i*4], D[0,t*4]] in Pair): 
                       Pair[k,:] =  [D[0,i*4], D[0,t*4]] 
                       k = k+1
    return Pair

Data_Mediterrean = np.genfromtxt('Data/MedSeaIDs.txt', delimiter=',')

Pair = SelectPair(Data_Mediterrean) #Find chance pairs

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

np.savetxt('Data/Pair.txt',Pair,delimiter=',')  #Save chance pair IDs
np.savetxt('Data/UnPair.txt',UnPair,delimiter=',')

# -*- coding: utf-8 -*-
"""
Chance pair distance timeseries

Created on Sat Oct 12 13:42:52 2019

@author: Gebruiker
"""
import numpy as np
import pandas as pd

def ComputeDistance(ID1,ID2,Data_Mediterrenean):
    id1 = []
    id2 = []
    for i in range(len(Data_Mediterrenean[0])):
        if Data_Mediterrenean[0,i] == ID1:
            id1 +=[[Data_Mediterrenean[1,i],Data_Mediterrenean[2,i],Data_Mediterrenean[3,i]]]
        if Data_Mediterrenean[0,i] == ID2:
            id2 +=[[Data_Mediterrenean[1,i],Data_Mediterrenean[2,i],Data_Mediterrenean[3,i]]]
    id1 = np.asarray(id1)
    id2 = np.asarray(id2)
    distance = []
    time = []
    for i in range(len(id1)):
        for j in range(len(id2)):
            if id1[i,2]==id2[j,2]:
                distance += [np.sqrt((id1[i,0]-id2[j,0])**2+(id1[i,1]-id2[j,1])**2)]
                time += [id1[i,2]]
    return distance,time

if __name__ == "__main__":
   # d,t = ComputeDistance(62835700,62835790,Data_Mediterrean)
    pairs = np.genfromtxt('UnPair.txt', delimiter=',')
    for i in range(len(pairs)):
        d,t = ComputeDistance(pairs[i,0],pairs[i,1],Data_Mediterrean)
        np.savetxt('PairDistances/Pair{0}.csv'.format(i),np.asarray((d,t)),delimiter = ',')
        
    

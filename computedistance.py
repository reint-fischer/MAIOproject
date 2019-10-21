# -*- coding: utf-8 -*-
"""
Chance pair distance timeseries

Created on Sat Oct 12 13:42:52 2019

@author: Gebruiker
"""
import numpy as np
import pandas as pd

def ComputeDistance(ID1,ID2,Data_Mediterrenean):
    id1 = [] #select only the 1st ID from all Mediterrenean data
    id2 = [] #select only the 1st ID from all Mediterrenean data
    for i in range(len(Data_Mediterrenean[0])):
        if Data_Mediterrenean[0,i] == ID1: #select right ID
            id1 +=[[Data_Mediterrenean[1,i],Data_Mediterrenean[2,i],Data_Mediterrenean[3,i]]] #save latitude, longitude, time
        if Data_Mediterrenean[0,i] == ID2: #select right ID
            id2 +=[[Data_Mediterrenean[1,i],Data_Mediterrenean[2,i],Data_Mediterrenean[3,i]]] #save latitude, longitude, time
    id1 = np.asarray(id1) #save as array for easy indexing
    id2 = np.asarray(id2) #save as array for easy indexing
    distance = [] #generate empty distance timeseries
    time = [] #generate corresponding timeaxis
    for i in range(len(id1)): #compare all measurement data
        for j in range(len(id2)):
            if id1[i,2]==id2[j,2]: # if the time is equal
                distance += [np.sqrt((id1[i,0]-id2[j,0])**2+(id1[i,1]-id2[j,1])**2)*1000] #compute distance in meters and add to timeseries
                time += [id1[i,2]] #add timestamp to timeaxis
    mind = distance.index(min(distance)) #find the index of the minimum separation distance to slice both 'distance' and 'time'
    d1 = list(reversed(distance[:mind+1])) #slice the timeseries up to the minimum and reverse it to create a backward timeseries
    d2 = distance[mind:] #slice the timeseries from the minimum onwards to create a forward timeseries
    t1 = list(reversed(time[:mind+1])) #slice te timeaxis in the same way as the timeseries
    t2 = time[mind:] #slice te timeaxis in the same way as the timeseries
    for n in range(len(t1)-1): #check for continuity
        if t1[n]-1 != t1[n+1]: #In backward timeaxis each next timestep should be 1 smaller
            t1 = t1[:n] #slice continuous timeaxis
            d1 = d1[:n] #slice corresponding backward distance timeseries
            break #stop for-loop when discontinuity is found
    for n in range(len(t2)-1): #do the same for the forward timeseries
        if t2[n]+1 != t2[n+1]:
            t2 = t2[:n]
            d2 = d2[:n]
            break
    return distance,time,d1,d2,t1,t2,mind

if __name__ == "__main__":
    nd = np.genfromtxt('Data/MedSeaIDs.txt',delimiter=',')
#    d,t,d1,d2,t1,t2,mind = ComputeDistance(131969,131970,nd)
    pairs = np.genfromtxt('Data/UnPair.txt', delimiter=',')
    for i in range(len(pairs)):
        d,t,d1,d2,t1,t2,mind = ComputeDistance(pairs[i,0],pairs[i,1],nd)
        np.savetxt('Data/BackwardsDistances/BDPair{0}.csv'.format(i),np.asarray((d1,t1)),delimiter = ',')
        np.savetxt('Data/ForwardDistances/FDPair{0}.csv'.format(i),np.asarray((d2,t2)),delimiter = ',')  
#    

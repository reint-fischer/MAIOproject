# -*- coding: utf-8 -*-
"""
Plot chance pairs

Created on Fri Oct 11 14:46:23 2019

@author: Gebruiker
"""
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt

nPair = 32

d,t = np.genfromtxt('PairDistances/Pair{0}.csv'.format(nPair),delimiter = ',')
pairs = np.genfromtxt('Data/UnPair.txt', delimiter=',')
nd = np.genfromtxt('Data/MedSeaIDslonlat.txt', delimiter=',')
#%%
ID1 = []
for i in range(len(nd[0])):
    if nd[0,i] == pairs[nPair,0]:
        ID1 +=[nd[1,i],nd[2,i],nd[3,i]]
ID1 = np.array(ID1)
ID1 = np.reshape(ID1,(int(len(ID1)/3),3))

ID2 = []
for i in range(len(nd[0])):
    if nd[0,i] == pairs[nPair,1]:
        ID2 +=[nd[1,i],nd[2,i],nd[3,i]]
ID2 = np.array(ID2)
ID2 = np.reshape(ID2,(int(len(ID2)/3),3))

f1 = plt.figure(1,figsize=(9,4))
ax1 = plt.axes(projection=ccrs.PlateCarree())
ax1.add_feature(cfeature.NaturalEarthFeature('physical','land','50m',edgecolor='face',facecolor='k'))
ax1.set_extent([-5,36,30,45],crs=ccrs.PlateCarree())
p1 = ax1.scatter(ID1[:,1],ID1[:,0])
p2 = ax1.scatter(ID2[:,1],ID2[:,0])

f2 = plt.figure(2,figsize=(9,4))
ax2 = plt.axes()
ax2.plot(t,d)


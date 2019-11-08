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
import cartopy.mpl.ticker as cticker
import matplotlib.animation as anim
import seaborn as sns
import random
import matplotlib

pairs = np.genfromtxt('Data/UnPair.txt', delimiter=',')
nd = np.genfromtxt('Data/MedSeaIDslonlat.txt', delimiter=',')


#----- Plot trajectories and timeseries of a single pair
nPair = 46

b,tb = np.genfromtxt('Data/BackwardsDistances/BDPair{0}.csv'.format(nPair),delimiter = ',')
f,tf = np.genfromtxt('Data/ForwardDistances/FDPair{0}.csv'.format(nPair),delimiter = ',')
ID1f = []
for i in range(len(nd[0])):
    if nd[0,i] == pairs[nPair,0]:
        if nd[3,i] in tf:
            ID1f +=[nd[1,i],nd[2,i],nd[3,i]]
ID1f = np.array(ID1f)
ID1f = np.reshape(ID1f,(int(len(ID1f)/3),3))

ID2f = []
for i in range(len(nd[0])):
    if nd[0,i] == pairs[nPair,1]:
        if nd[3,i] in tf:
            ID2f +=[nd[1,i],nd[2,i],nd[3,i]]
ID2f = np.array(ID2f)
ID2f = np.reshape(ID2f,(int(len(ID2f)/3),3))

ID1b = []
for i in range(len(nd[0])):
    if nd[0,i] == pairs[nPair,0]:
        if tb.size>1:
            if nd[3,i] in tb:
                ID1b +=[nd[1,i],nd[2,i],nd[3,i]]
        else:
            if nd[3,i] in [tb]:
                ID1b +=[nd[1,i],nd[2,i],nd[3,i]]
ID1b = np.array(ID1b)
ID1b = np.reshape(ID1b,(int(len(ID1b)/3),3))

ID2b = []
for i in range(len(nd[0])):
    if nd[0,i] == pairs[nPair,1]:
        if tb.size>1:
            if nd[3,i] in tb:
                ID2b +=[nd[1,i],nd[2,i],nd[3,i]]
        else:
            if nd[3,i] in [tb]:
                ID2b +=[nd[1,i],nd[2,i],nd[3,i]]
ID2b = np.array(ID2b)
ID2b = np.reshape(ID2b,(int(len(ID2b)/3),3))

f1 = plt.figure(1,figsize=(9,9))
ax1 = plt.axes(projection=ccrs.PlateCarree())
ax1.add_feature(cfeature.NaturalEarthFeature('physical','land','50m',edgecolor='face',facecolor='k'))
ax1.add_feature(cfeature.NaturalEarthFeature('physical','ocean','50m',edgecolor='face',facecolor='c'))
ax1.set_extent([-5,36,30,45],crs=ccrs.PlateCarree())
ax1.set_xticks([0,10,20,30], crs=ccrs.PlateCarree())
ax1.set_xticklabels([0,10,20,30])
ax1.set_yticks([30,35,40,45], crs=ccrs.PlateCarree())
ax1.set_yticklabels([30,35,40,45])
lon_formatter = cticker.LongitudeFormatter()
lat_formatter = cticker.LatitudeFormatter()
ax1.xaxis.set_major_formatter(lon_formatter)
ax1.yaxis.set_major_formatter(lat_formatter)
def animate(i):
    ax1.plot(ID1f[:i,1],ID1f[:i,0],transform=ccrs.Geodetic(),color='red')
    ax1.plot(ID2f[:i,1],ID2f[:i,0],transform=ccrs.Geodetic(),color='red')
    ax1.plot(list(reversed(ID1b[:,1]))[:i],list(reversed(ID1b[:,0]))[:i],transform=ccrs.Geodetic(),color='blue')
    ax1.plot(list(reversed(ID2b[:,1]))[:i],list(reversed(ID2b[:,0]))[:i],transform=ccrs.Geodetic(),color='blue')
ani = anim.FuncAnimation(f1, animate)
#ani.save('Figures/anim2.mp4', fps=3, extra_args=['-vcodec', 'libx264'])

f4 = plt.figure(4,figsize=(9,4))
ax4 = plt.axes()
ax4.plot(b)

f5 = plt.figure(5,figsize=(9,4))
ax5 = plt.axes()
ax5.plot(f)

# Plot all trajectories
palette = sns.color_palette("hls", len(pairs))
f7 = plt.figure(7,figsize=(9,4))
ax7 = plt.axes(projection=ccrs.PlateCarree())
ax7.add_feature(cfeature.NaturalEarthFeature('physical','land','50m',edgecolor='face',facecolor='k'))
ax7.set_extent([-5,36,30,45],crs=ccrs.PlateCarree())
ax7.set_xticks([-5,0,5,10,15,20,25,30,35], crs=ccrs.PlateCarree())
ax7.set_xticklabels([-5,0,5,10,15,20,25,30,35])
ax7.set_yticks([30,35,40,45], crs=ccrs.PlateCarree())
ax7.set_yticklabels([30,35,40,45])
lon_formatter = cticker.LongitudeFormatter()
lat_formatter = cticker.LatitudeFormatter()
ax7.xaxis.set_major_formatter(lon_formatter)
ax7.yaxis.set_major_formatter(lat_formatter)
for h in range(len(pairs)):
    ID1 = []
    for i in range(len(nd[0])):
        if nd[0,i] == pairs[h,0]:
            ID1 +=[nd[1,i],nd[2,i],nd[3,i]]
    ID1 = np.array(ID1)
    ID1 = np.reshape(ID1,(int(len(ID1)/3),3))
    
    ID2 = []
    for i in range(len(nd[0])):
        if nd[0,i] == pairs[h,1]:
            ID2 +=[nd[1,i],nd[2,i],nd[3,i]]
    ID2 = np.array(ID2)
    ID2 = np.reshape(ID2,(int(len(ID2)/3),3))
    p1 = ax7.plot(ID1[:,1],ID1[:,0],color = palette[h],transform=ccrs.Geodetic())
    p2 = ax7.plot(ID2[:,1],ID2[:,0],color = palette[h],transform=ccrs.Geodetic())

##%% Plot all trajectories with different styles for forward and backward
palette = sns.hls_palette(len(pairs), l=.7, s=.8)
palette2 = sns.hls_palette(len(pairs), l=.3, s=.8)
f8 = plt.figure(8,figsize=(9,4))
f8.tight_layout()
ax8 = plt.axes(projection=ccrs.PlateCarree())
ax8.add_feature(cfeature.NaturalEarthFeature('physical','land','50m',edgecolor='face',facecolor='k'))
ax8.set_extent([-5,36,30,45],crs=ccrs.PlateCarree())
ax8.set_xticks([-5,0,5,10,15,20,25,30,35], crs=ccrs.PlateCarree())
ax8.set_xticklabels([-5,0,5,10,15,20,25,30,35])
ax8.set_yticks([30,35,40,45], crs=ccrs.PlateCarree())
ax8.set_yticklabels([30,35,40,45])
lon_formatter = cticker.LongitudeFormatter()
lat_formatter = cticker.LatitudeFormatter()
ax8.xaxis.set_major_formatter(lon_formatter)
ax8.yaxis.set_major_formatter(lat_formatter)
for h in range(len(pairs)):
    b,tb = np.genfromtxt('Data/BackwardsDistances/BDPair{0}.csv'.format(h),delimiter = ',')
    f,tf = np.genfromtxt('Data/ForwardDistances/FDPair{0}.csv'.format(h),delimiter = ',')
    ID1f = []
    for i in range(len(nd[0])):
        if nd[0,i] == pairs[h,0]:
            if nd[3,i] in tf:
                ID1f +=[nd[1,i],nd[2,i],nd[3,i]]
    ID1f = np.array(ID1f)
    ID1f = np.reshape(ID1f,(int(len(ID1f)/3),3))
    
    ID2f = []
    for i in range(len(nd[0])):
        if nd[0,i] == pairs[h,1]:
            if nd[3,i] in tf:
                ID2f +=[nd[1,i],nd[2,i],nd[3,i]]
    ID2f = np.array(ID2f)
    ID2f = np.reshape(ID2f,(int(len(ID2f)/3),3))
    
    ID1b = []
    for i in range(len(nd[0])):
        if nd[0,i] == pairs[h,0]:
            if tb.size>1:
                if nd[3,i] in tb:
                    ID1b +=[nd[1,i],nd[2,i],nd[3,i]]
            else:
                if nd[3,i] in [tb]:
                    ID1b +=[nd[1,i],nd[2,i],nd[3,i]]
    ID1b = np.array(ID1b)
    ID1b = np.reshape(ID1b,(int(len(ID1b)/3),3))
    
    ID2b = []
    for i in range(len(nd[0])):
        if nd[0,i] == pairs[h,1]:
            if tb.size>1:
                if nd[3,i] in tb:
                    ID2b +=[nd[1,i],nd[2,i],nd[3,i]]
            else:
                if nd[3,i] in [tb]:
                    ID2b +=[nd[1,i],nd[2,i],nd[3,i]]
    ID2b = np.array(ID2b)
    ID2b = np.reshape(ID2b,(int(len(ID2b)/3),3))

    p1 = ax8.plot(ID1f[:,1],ID1f[:,0],transform=ccrs.Geodetic(),color = palette[h])
    p2 = ax8.plot(ID2f[:,1],ID2f[:,0],transform=ccrs.Geodetic(),color = palette[h])
    p3 = ax8.plot(ID1b[:,1],ID1b[:,0],transform=ccrs.Geodetic(),color = palette2[h],linestyle=':')
    p4 = ax8.plot(ID2b[:,1],ID2b[:,0],transform=ccrs.Geodetic(),color = palette2[h],linestyle=':')

f9 = plt.figure(9,figsize=(9,4))
f9.tight_layout()
ax9 = plt.axes(projection=ccrs.PlateCarree())
ax9.add_feature(cfeature.NaturalEarthFeature('physical','land','50m',edgecolor='face',facecolor='k'))
ax9.set_extent([-5,36,30,45],crs=ccrs.PlateCarree())
ax9.set_xticks([-5,0,5,10,15,20,25,30,35], crs=ccrs.PlateCarree())
ax9.set_xticklabels([-5,0,5,10,15,20,25,30,35])
ax9.set_yticks([30,35,40,45], crs=ccrs.PlateCarree())
ax9.set_yticklabels([30,35,40,45])
lon_formatter = cticker.LongitudeFormatter()
lat_formatter = cticker.LatitudeFormatter()
ax9.xaxis.set_major_formatter(lon_formatter)
ax9.yaxis.set_major_formatter(lat_formatter)
ax9.scatter(nd[2],nd[1],s=0.1)
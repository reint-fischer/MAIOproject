# -*- coding: utf-8 -*-
"""
Opening netCDF drifter file and saving Mediterranean data

Created on Wed Oct  9 10:22:58 2019

@author: Gebruiker
"""
from netCDF4 import Dataset,num2date
import numpy as np
import xarray as xr

def OpenDrifterdata(filelocation): #Open netcdf file and unpack relevant parameters
    file = Dataset(filelocation,'r+',format = 'NETCDF4')
    time = file.variables['TIME'][:] #timestamp of GPS location
    latitude = file.variables['LAT'][:] #latitude of corresponding drifter ID at corresponding time
    longitude = file.variables['LON'][:] #longitude of corresponding drifter ID at corresponding time
    ID = file.variables['ID'][:] 
    file.close()
    return ID,latitude,longitude,time

def selectMedsea(ID,lat,lon,time):
    newID = []
    newlat = []
    newlon = []
    newtime = []
    for i in range(len(lat)):
        if lat[i]<45 and lat[i]>30 and lon[i]<36 and lon[i]>5:
            newID += [ID[i]]
            newlat += [lat[i]]
            newlon += [lon[i]]
            newtime += [time[i]]
    return np.array([newID,newtime,newlat,newlon])


#%%
if __name__ == '__main__':
    ID,latitude,longitude,time = OpenDrifterdata('C:/Users/Gebruiker/Downloads/driftertrajGPS_1.03.nc')
    #%%
    nd = selectMedsea(ID,latitude,longitude,time)
#    np.savetxt('MedSeaIDs.csv',nd,delimiter=',')

from netCDF4 import Dataset,num2date
import numpy as np
import math as m
import pandas as pd

#-----------------------FUNCTIONS--------------------

#---> Convert the data set of lat and lon to distance from the left down angle (30N:-5E)
def Conversion(D):  
    Earth_radius = 6371.0   # [km]
    for i in range(len(D[1,:])) :
        if D[2,i]>36:                 # Used for having data in the strait of Gilbilter where LON = 355
            D[2,i]=D[2,i]-360         # In this way it has a negative value for instance 355 -> -5E and so the new D[2,i] = 0  
        D[2,i] = 2*m.pi*Earth_radius*m.cos(np.deg2rad(D[1,i]))*(D[2,i]+5)/360
        D[1,i] = 2*m.pi*Earth_radius*(D[1,i] - 30)/360 
    return D    
      
#---> Open netcdf file and unpack relevant parameters
def OpenDrifterdata(filelocation): 
    file = Dataset(filelocation,'r+',format = 'NETCDF4')
    time = file.variables['TIME'][:]          # timestamp of GPS location
    latitude = file.variables['LAT'][:]       # latitude of corresponding drifter ID at corresponding time
    longitude = file.variables['LON'][:]      # longitude of corresponding drifter ID at corresponding time
    ID = file.variables['ID'][:] 
    return ID,latitude,longitude,time

#---> Selecting the mediterrean data
def selectMedsea(ID,lat,lon,time):
    newID = []
    newlat = []
    newlon = []
    newtime = []
    for t in range(len(lat)): 
        if lat[t]<45 and lat[t]>30 and lon[t]<36: # for lon bigger that the London longitude
            newID += [ID[t]]
            newlat += [lat[t]]
            newlon += [lon[t]]
            newtime += [time[t]]
        if lat[t]<40 and lat[t]>33 and lon[t]>355: #for lon smaller that hte London longitude
            newID += [ID[t]]
            newlat += [lat[t]]
            newlon += [lon[t]]
            newtime += [time[t]]
    return np.array((newID,newlat,newlon,newtime))

#-----------------------MAIN--------------------

if __name__ == '__main__':
    ID,latitude,longitude,time = OpenDrifterdata('/home/giovanni/MAIO/main/driftertrajGPS_1.03.nc')    # Open netcdf file [FUNCTION 2]
    Data_Mediterrean = selectMedsea(ID,latitude, longitude, time)                                      # Select measurements [FUNCTION 3]
    Data_Mediterrean = Conversion(Data_Mediterrean)                                                    # Reanaysis distance [FUNCTION 1]
    np.savetxt('MedSeaIDs.txt',Data_Mediterrean,delimiter=',') 					       # Save data (reanalysis with distances at reference)	

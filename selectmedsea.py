from netCDF4 import Dataset,num2date
import numpy as np
import math as m
import pandas as pd

def Conversion(D):  #Convert the data set to distance from the left down angle (30N 5E)
    Earth_radius = 6371.0 #[m]
    for i in range(len(D[1,:])):
        if D[2,i]>36:
            D[2,i]=D[2,i]-360
        D[2,i] = 2*m.pi*Earth_radius*m.cos(np.deg2rad(D[1,i]))*(D[2,i]+5)/360
        D[1,i] = 2*m.pi*Earth_radius*(D[1,i] - 30)/360 
    return D    
      

def OpenDrifterdata(filelocation): #Open netcdf file and unpack relevant parameters
    file = Dataset(filelocation,'r+',format = 'NETCDF4')
    time = file.variables['TIME'][:] #timestamp of GPS location
    latitude = file.variables['LAT'][:] #latitude of corresponding drifter ID at corresponding time
    longitude = file.variables['LON'][:] #longitude of corresponding drifter ID at corresponding time
    ID = file.variables['ID'][:] 
    return ID,latitude,longitude,laterr,lonerr,time

def selectMedsea(ID,lat,lon,time):
    newID = []
    newlat = []
    newlon = []
    newtime = []
    for t in range(len(lat)): 
        if lat[t]<45 and lat[t]>30 and lon[t]<36: 
            newID += [ID[t]]
            newlat += [lat[t]]
            newlon += [lon[t]]
            newtime += [time[t]]
        if lat[t]<40 and lat[t]>33 and lon[t]>355: 
            newID += [ID[t]]
            newlat += [lat[t]]
            newlon += [lon[t]]
            newtime += [time[t]]
    return np.array((newID,newlat,newlon,newtime))
#%%
if __name__ == '__main__':
    ID,latitude,longitude,laterr,lonerr,time = OpenDrifterdata('C:/Users/Gebruiker/Downloads/driftertrajGPS_1.03.nc') #Open netcdf file
#    Data_Mediterrean = selectMedsea(ID,latitude, longitude, time) #select measurements made in Mediterrenean
#    np.savetxt('Data/MedSeaIDslonlat.txt',Data_Mediterrean,delimiter=',') #save data of all measurements made in the Mediterrenean
#    Data_Mediterrean = Conversion(Data_Mediterrean) #convert measurements to flat grid
#    np.savetxt('Data/MedSeaIDs.txt',Data_Mediterrean,delimiter=',') #save data of all converted measurements

    

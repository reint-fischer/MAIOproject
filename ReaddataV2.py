from netCDF4 import Dataset,num2date
import numpy as np
import math as m

def Conversion(D):  #Convert the data set to distance from the left down angle (30N 5E)
    Earth_radius = 6371.0 #[m]
    for i in range(len(D[1,:])): 
       D[2,i] = 2*m.pi*Earth_radius*m.cos(np.deg2rad(D[1,i]))*(D[2,i]-5)/360
       D[1,i] = 2*m.pi*Earth_radius*(D[1,i] - 30)/360 
    return D    

def SelectPair(D):
    Pair=[]
    for i in range ( int(len(D[2,:])/4) ): 
    print( 100*(i*1.0)/(1.0*int(len(D[2,:])/4)) '%' )
        for t in range ( int(len(D[2,:])/4) ):
            if D[0,i*4]!= D[0,t*4] : 
                space_distance = (D[2,i*4]-D[2,t*4])*(D[2,i*4]-D[2,t*4]) + (D[1,i*4]-D[1,t*4])*(D[1,i*4]-D[1,t*4])
                time_distance = (D[3,i*4] - D[3,t*4])*(D[3,i*4] - D[3,t*4]) 
                if space_distance < 100 and time_distance < 36 :
                   if not( [D[0,i*4], D[0,t*4]] in Pair): Pair +=  [D[0,i*4], D[0,t*4]] 
    return Pair 

def OpenDrifterdata(filelocation): #Open netcdf file and unpack relevant parameters
    file = Dataset(filelocation,'r+',format = 'NETCDF4')
    time = file.variables['TIME'][:] #timestamp of GPS location
    latitude = file.variables['LAT'][:] #latitude of corresponding drifter ID at corresponding time
    longitude = file.variables['LON'][:] #longitude of corresponding drifter ID at corresponding time
    ID = file.variables['ID'][:] 
    return ID,latitude,longitude,time

def selectMedsea(ID,lat,lon,time):
    newID = []
    newlat = []
    newlon = []
    newtime = []
    for t in range(len(lat)): 
        if lat[t]<45 and lat[t]>30 and lon[t]<36 and lon[t]>-5 : 
            newID += [ID[t]]
            newlat += [lat[t]]
            newlon += [lon[t]]
            newtime += [time[t]] 
    return np.array((newID,newlat,newlon,newtime))

if __name__ == '__main__':
    ID,latitude,longitude,time = OpenDrifterdata('/home/giovanni/MAIO/Python/driftertrajGPS_1.03.nc')
    Data_Mediterrean = selectMedsea(ID,latitude, longitude, time)
    Data_Mediterrean = Conversion(Data_Mediterrean)
    Pair = SelectPair(Data_Mediterrean)
    np.savetxt('MedSeaIDs.txt',Data_Mediterrean,delimiter=',')
    np.savetxt('Pair.txt',Pair,delimiter=',') 
    

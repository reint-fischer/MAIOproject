# -*- coding: utf-8 -*-
"""
Opening netCDF drifter file and saving Mediterranean data

Created on Wed Oct  9 10:22:58 2019

@author: Gebruiker
"""
from netCDF4 import Dataset,num2date

def OpenDrifterdata(filelocation):
    file = Dataset(filelocation,'r+',format = 'NETCDF4')
    time = num2date(file.variables['time'][:],file.variables['time'].units)
    latitude = file.variables['lat'][:]
    longitude = file.variables['lon'][:]
    ID = file.variables['ID'][:]


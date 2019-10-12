# -*- coding: utf-8 -*-
"""
Opening netCDF drifter file and saving Mediterranean data using XARRAY

Created on Thu Oct 10 14:30:28 2019

@author: reint fischer
"""

import numpy as np
import xarray as xr
import pandas as pd

def reshapedrifterdata(filestring):
    DS = xr.open_dataset(filestring)
    DS = DS.drop('U')
    DS = DS.drop('V')
    DS = DS.drop('LAT_ERR')
    DS = DS.drop('LON_ERR')
    DS = DS.drop('U_ERR')
    DS = DS.drop('V_ERR')
    DS = DS.drop('GAP')
    DS = DS.drop('DROGUE')
    DS= DS.set_coords('ID')
    newID = []
    newtime = []
    for i in range(len(DS.ID.values)):
        if not(DS.ID.values[i] in newID):
            newID += [DS.ID.values[i]]
        if not(DS.TIME.values[i] in newtime):
            newtime += [DS.TIME.values[i]]
    return DS
#%%
if __name__ == "__main__":
    DS = xr.open_dataset('C:/Users/Gebruiker/Downloads/driftertrajGPS_1.03.nc')
    DS = DS.drop('U')
    DS = DS.drop('V')
    DS = DS.drop('LAT_ERR')
    DS = DS.drop('LON_ERR')
    DS = DS.drop('U_ERR')
    DS = DS.drop('V_ERR')
    DS = DS.drop('GAP')
    DS = DS.drop('DROGUE')
    DS= DS.set_coords('ID')
    x = DS.TIME.values
    y = DS.ID.values
#    ind = pd.MultiIndex.from_product((x,y),names=('TIME','ID'))
    
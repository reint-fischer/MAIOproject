# -*- coding: utf-8 -*-
"""
Opening netCDF drifter file and saving Mediterranean data using XARRAY

Created on Thu Oct 10 14:30:28 2019

@author: reint fischer
"""

import numpy as np
import xarray as xr
import pandas as pd

def SelectMedSea(filestring):
    DS = xr.open_dataset(filestring)
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
    ind = pd.MultiIndex.from_product((x,y),names=('TIME','ID'))
    
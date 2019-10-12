# -*- coding: utf-8 -*-
"""
Opening netCDF drifter file and saving Mediterranean data using XARRAY

Created on Thu Oct 10 14:30:28 2019

@author: reint fischer
"""

import numpy as np
import xarray as xr
import pandas as pd

def reshapedrifterdata(meddata):
    ID = meddata[0]
    time = meddata[1]
    lat = meddata[2]
    lon = meddata[3]
    
    newID = []
    newtime = []
    for i in range(len(ID)):
        if not(ID[i] in newID):
            newID += [ID[i]]
        if not(time[i] in newtime):
            newtime += [time[i]]
    nd = np.full((len(newtime),len(newID),2),np.nan)
    for j in range(len(lat)):
        for u in range(len(newtime)):
            for v in range(len(newID)):
                if ID[j]==newID[v] and time[j]==newtime[u]:
                    nd[u,v,0] = lat[j]
                    nd[u,v,1] = lon[j]
    return nd
#%%
if __name__ == "__main__":
    redata = reshapedrifterdata(nd)
#    ind = pd.MultiIndex.from_product((x,y),names=('TIME','ID'))
    
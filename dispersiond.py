# -*- coding: utf-8 -*-
"""
D2 Ensemble

Created on Wed Oct 16 22:26:40 2019

@author: Gebruiker
"""
import numpy as np
import matplotlib.pyplot as plt

pairs = np.genfromtxt('Data/UnPair.txt', delimiter=',')
b = []
f = []
f1 = plt.figure(1)
ax1 = plt.axes()
#plt.yscale('log')
f2 = plt.figure(2)
ax2 = plt.axes()
#plt.yscale('log')
for i in range(len(pairs)):
    b += [np.square(np.genfromtxt('BackwardsDistances/BDPair{0}.csv'.format(i),delimiter = ',')[0])]
    f += [np.square(np.genfromtxt('ForwardDistances/FDPair{0}.csv'.format(i),delimiter = ',')[0])]
    ax1.loglog(b[i])
    ax2.loglog(f[i])



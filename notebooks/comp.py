# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 17:01:25 2016

@author: schiavon
"""

from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt

t = loadmat('../data/psi_comp.mat')
lsm = t['ls']
lim = t['li']
psim = t['psi']

a = np.load('../data/psi_comp.npz')
ls = a['ls']
li = a['li']
psi = a['psi']

plt.figure()
plt.imshow((np.abs(psim)-psi[0].T)/np.abs(psim),interpolation='nearest')

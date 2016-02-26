# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 11:01:00 2016

@author: schiavon
"""

from __future__ import print_function,division

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import bennink as B
import myintegrate
from scipy.optimize import minimize,fmin

from joblib import Parallel,delayed
import time
    
Lambda = 10e-6
L = 0.03

axp = axi = 'y'
axs = 'z'

ls = li = 810e-9
lp = 405e-9

wp = np.arange(10,80,1)*1e-6

f = lambda a,lp,ls,li,wp,axp,axs,axi,L,Lambda: -np.abs(B.O(lp,ls,li,wp,a*wp,a*wp,axp,axs,axi,L,Lambda))

#def funcycle(wp):
#    ret = 
#    return ret[0]
#
#r = Parallel(n_jobs=8)(delayed(funcycle)(i) for i in range(wp.size))

psi = np.zeros(wp.size)
a = np.zeros(wp.size)

for i in range(wp.size):
    t1 = time.time()
    a[i],psi[i],c,d,e = fmin(f,2,args=(lp,ls,li,wp[i],axp,axs,axi,L,Lambda),full_output=1)
    t2 = time.time()
    print('Cycle',i,'in',t2-t1,'s')
    
np.savez('../data/psi_max_wp.npz',lp=lp,ls=ls,li=li,axp=axp,axs=axp,axi=axi,Lambda=Lambda,L=L,wp=wp,psi=psi,a=a)

#%% plot results
plt.figure()
plt.plot(wp*1e6,-psi)
plt.xlabel('wp [um]')
plt.ylabel(r'Bennink $\psi$ [a.u.]')
plt.title('Spatial overlap')


plt.figure()
plt.plot(wp*1e6,a*wp*1e6)
plt.xlabel('wp [um]')
plt.ylabel('ws,wi [um]')
plt.title('Waist maximizing the spatial overlap')
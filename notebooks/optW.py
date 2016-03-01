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
#import myintegrate
from scipy.optimize import fmin
from scipy.constants import c

from mpl_toolkits.mplot3d import Axes3D

import os.path

#from joblib import Parallel,delayed
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

if not os.path.isfile('../data/psi_max_wp.npz'):
    for i in range(wp.size):
        t1 = time.time()
        a[i],psi[i],c,d,e = fmin(f,2,args=(lp,ls,li,wp[i],axp,axs,axi,L,Lambda),full_output=1)
        t2 = time.time()
        print('Cycle',i,'in',t2-t1,'s')
        
    np.savez('../data/psi_max_wp.npz',lp=lp,ls=ls,li=li,axp=axp,axs=axp,axi=axi,Lambda=Lambda,L=L,wp=wp,psi=psi,a=a)
else:
    data = np.load('../data/psi_max_wp.npz')
    psi = data['psi']
    a = data['a']

#%% plot results
plt.close('all')

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

pos_max = np.argmin(psi)
print('Maximum psi(wp,ws) for wp =',wp[pos_max]*1e6,'um, with ws = wi =',wp[pos_max]*a[pos_max]*1e6,'um')
print('psi(omega_p/2,omega_p/2) =',-psi[pos_max])

#%% compute the psi function for the optimized parameters
wpm = wp[pos_max]
wsm = wim = wp[pos_max]*a[pos_max]

omega_p = 2*np.pi*c/lp
ls = np.arange(800,820,1)*1e-9
omega_s = 2*np.pi*c/ls
omega_i = np.ones(len(omega_s))*omega_p-omega_s
li = 2*np.pi*c/omega_i

O_omega = np.zeros( (omega_s.size,omega_i.size) )

for i in range(ls.size):
    for j in range(li.size):
        O_omega[i,j] = np.abs(B.O(lp,ls[i],li[j],wpm,wsm,wim,axp,axs,axi,L,Lambda))

#%%
fig = plt.figure()
ax = Axes3D(fig)

ls2,li2 = np.meshgrid(ls,li)
surf = ax.plot_surface(ls2,li2,O_omega,rstride=1, cstride=1, cmap=cm.coolwarm,linewidth=0, antialiased=False)

fig.colorbar(surf, shrink=0.5, aspect=5)

# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 11:01:00 2016

@author: schiavon
"""

from __future__ import print_function,division

import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.cm as cm
import bennink as B
#import myintegrate
from scipy.optimize import fmin
from scipy.constants import c

import os.path

from joblib import Parallel,delayed
import time
    
Lambda = 10e-6
L = 0.03

axp = axi = 'y'
axs = 'z'

ls = li = 809e-9
lp = 404.9e-9

wp = np.arange(10,80,1)*1e-6

f = lambda a,lp,ls,li,wp,axp,axs,axi,L,Lambda: -np.abs(B.O(lp,ls,li,wp,a*wp,a*wp,axp,axs,axi,L,Lambda))

#def funcycle(wp):
#    ret = 
#    return ret[0]
#
#r = Parallel(n_jobs=8)(delayed(funcycle)(i) for i in range(wp.size))

O = np.zeros(wp.size)
a = np.zeros(wp.size)

if not os.path.isfile('../data/psi_max_wp.npz'):
    pass
    t1 = time.time()
    def func(i):    
        a,o,c,d,e = fmin(f,2,args=(lp,ls,li,wp[i],axp,axs,axi,L,Lambda),full_output=1)
        return a,o

    a,O = Parallel(n_jobs=8)(delayed(func)(i) for i in range(wp.size))

    t2 = time.time()

        
    np.savez('../data/psi_max_wp.npz',lp=lp,ls=ls,li=li,axp=axp,axs=axp,axi=axi,Lambda=Lambda,L=L,wp=wp,O=O,a=a)
else:
    data = np.load('../data/psi_max_wp.npz')
    O = data['O']
    a = data['a']

#%% plot results
plt.close('all')

plt.figure()
plt.plot(wp*1e6,-O)
plt.xlabel('wp [um]')
plt.ylabel(r'Bennink $\psi$ [a.u.]')
plt.title('Spatial overlap')


plt.figure()
plt.plot(wp*1e6,a*wp*1e6)
plt.xlabel('wp [um]')
plt.ylabel('ws,wi [um]')
plt.title('Waist maximizing the spatial overlap')

pos_max = np.argmin(O)
print('Maximum psi(wp,ws) for wp =',wp[pos_max]*1e6,'um, with ws = wi =',wp[pos_max]*a[pos_max]*1e6,'um')
print('psi(omega_p/2,omega_p/2) =',-O[pos_max])

#%% compute the psi function for the optimized parameters
wpm = wp[pos_max]
wsm = wim = wp[pos_max]*a[pos_max]

omega_p = 2*np.pi*c/lp
lsc = 809.0e-9
steps = 10
stepsize = 1e-11

ls = np.linspace(lsc - steps/2*stepsize, ls + steps/2*stepsize, steps)

lic = 1/(1/lp - 1/lsc)
li = np.linspace(lic - steps/2*stepsize, lic + steps/2*stepsize, steps)

O_omega = np.zeros( (ls.size,li.size) )

t1 = time.time()
def int_cycle(i):
   ret = np.zeros(li.size)
   for j in range(li.size):
       ret[j] = np.abs(B.O(lp,ls[i],li[j],wpm,wsm,wim,axp,axs,axi,L,Lambda))
   return ret

r = Parallel(n_jobs=8)(delayed(int_cycle)(i) for i in range(ls.size))

for i in range(ls.size):
    O_omega[:,i] = r[i]
    
t2 = time.time()

print('O_omega in the interval ('+str(ls[0])+','+str(ls[-1])+') computed in',t2-t1,'s')

#%%
plt.figure()

plt.imshow(O_omega,origin='lower',extent=(ls[0]*1e9,ls[-1]*1e9,li[0]*1e9,li[-1]*1e9),aspect='auto',vmin=0,interpolation='none')
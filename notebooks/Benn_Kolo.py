# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import bennink as B
import myintegrate

from joblib import Parallel,delayed
import time
    
Lambda = 10e-6
wp = 35e-6
ws = wi = 50e-6
L = 0.03

axp = axi = 'y'
axs = 'z'

#%% Phasematching
li = ls = np.arange(809.9e-9,810.1e-9,5e-12)
#psi1 = np.zeros((ls.size,li.size))
Ns = np.array([500,1000,2000,10000])

psi = np.zeros((Ns.size,ls.size,li.size))

#t0 = time.time()
#for i in range(ls.size):
#    for j in range(li.size):
#        psi[j,i] = np.abs(psi_Bennink(1/(1/ls[i]+1/li[j]),ls[i],li[j],wp,ws,wi,axp,axs,axi,'stupid'))
#t1 = time.time()

for j in range(Ns.size):
    myintegrate.N = Ns[j]    
    
    t2 = time.time()
    def int_cycle(i):
        loc = np.zeros(li.size)
        for j in range(li.size):
            loc[j] = np.abs(B.psi(1/(1/ls[i]+1/li[j]),ls[i],li[j],wp,ws,wi,axp,axs,axi,L,Lambda,'stupid'))
        return loc
    
    r = Parallel(n_jobs=8)(delayed(int_cycle)(i) for i in range(ls.size))
    
    for i in range(ls.size):
        psi[j,:,i] = r[i]
    t3 = time.time()
    
    #print('Sequential in',t1-t0,'s')
    print('Parallel in',t3-t2,'s for N=',Ns[j])

    
#%% Figures
plt.figure()
plt.imshow(psi[3], origin='lower', extent=(ls.min()*1e9, ls.max()*1e9, li.min()*1e9, li.max()*1e9),interpolation='nearest', cmap=cm.gist_rainbow)
plt.colorbar()
plt.ticklabel_format(axis='both',style='plain',useOffset=False)
plt.grid()


#ax=plt.gca()
#fmt=matplotlib.ticker.ScalarFormatter(useOffset=False)
#fmt.set_scientific(False)
#ax.xaxis.set_major_formatter(fmt)
#ax.yaxis.set_major_formatter(fmt)

#%% Optimal grating period
#lp = 405e-9
#ls = li = 810e-9
#axp = axi = 'y'
#axs = 'z'
#
#Lopt = (n(lp,axp)/lp - n(ls,axs)/ls - n(li,axi)/li)**-1
#print('Grating period =',Lopt*1e6,'um')

#%%

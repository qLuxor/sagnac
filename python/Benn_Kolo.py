# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import scipy.integrate as integ
from sellmeier import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib

from joblib import Parallel,delayed
import time

c = 299792458 # [m/s]

#plt.close('all')



def n(l,axis):
    ''' Sellmeier equations for PPKTP (from Kato-Takaoka, at 20°C),
        x and y are ordinary axes, while z is the extraordinary one.
        l is the wavelength in m'''
    l = l*1e6 # transform wavelength to micron
    if axis == 'x':
        n = np.sqrt(3.291 + 0.0414/(l**2 - 0.03978) + 9.35522/(l**2 - 31.45571))
    elif axis == 'y':
        n = np.sqrt(3.45018 + 0.04341/(l**2 - 0.04597) + 16.98825/(l**2 - 39.43799))
    elif axis == 'z':
        n = np.sqrt(4.59423 + 0.06206/(l**2 - 0.04763) + 110.80672/(l**2 - 86.12171))
    return n
    
#def n(l,axis):
#    return sellmeier_ktp(l,axis)
    
Lambda = 10e-6
wp = 35e-6
ws = wi = 50e-6
L = 0.03

axp = axi = 'y'
axs = 'z'

def omega2lambda(omega):
    return 2*np.pi*c/omega

def lambda2omega(l):
    return 2*np.pi*c/l

omega_p = lambda2omega(405e-9)
os = lambda2omega(805e-9)
#oi = lambda2omega(810e-9)

def complex_integral(func,a,b,args,intype='stupid'):
    real_func = lambda z: np.real(func(z,*args))
    imag_func = lambda z: np.imag(func(z,*args))
    if intype == 'quad':
        real_int = integ.quad(real_func,a,b)
        imag_int = integ.quad(imag_func,a,b)
#        print(real_int)
#        print(imag_int)
        return real_int[0] + 1j * imag_int[0]
    elif intype == 'quadrature':
        real_int = integ.quadrature(real_func,a,b)
        imag_int = integ.quadrature(imag_func,a,b)
#        print(real_int)
#        print(imag_int)
        return real_int[0] + 1j * imag_int[0]
    elif intype == 'romberg':
        real_int = integ.romberg(real_func,a,b)
        imag_int = integ.romberg(imag_func,a,b)
#        print(real_int)
#        print(imag_int)
        return real_int + 1j * imag_int
    elif intype == 'stupid':
        Npoints = 500
        z,dz = np.linspace(a,b,Npoints,retstep=True)
        real_int = np.sum(real_func(z))*dz
        imag_int = np.sum(imag_func(z))*dz
#        print(real_int)
#        print(imag_int)
        return real_int + 1j*imag_int
        
#%% SPDC type-II -> processo y-z-y

# Bennink
def func_Bennink(z,lp,ls,li,wp,ws,wi,axp,axs,axi):
    omega_p = lambda2omega(lp)
    omega_s = lambda2omega(ls)
    omega_i = lambda2omega(li)
    
    def qp(z):
        return wp**2 + 1j*z*lp/np.pi/n(lp,axp)

    def qs(z):
        return ws**2 + 1j*z*ls/np.pi/n(ls,axs)
        
    def qi(z):
        return wi**2 + 1j*z*li/np.pi/n(li,axi)
    
    return wp*ws*wi*np.exp(1j*( (n(lp,axp)*omega_p - n(ls,axs)*omega_s - n(li,axi)*omega_i)/c + 2*np.pi/Lambda )*z) / (np.conj(qs(z))*np.conj(qi(z)) + qp(z)*np.conj(qi(z)) + qp(z)*np.conj(qs(z)))
#    return wp*ws*wi / (np.conj(qs(z))*np.conj(qi(z)) + qp(z)*np.conj(qi(z)) + qp(z)*np.conj(qs(z)))

def psi_Bennink(lp,ls,li,wp,ws,wi,axp,axs,axi,intype):
    return complex_integral(func_Bennink,-L/2,L/2,(lp,ls,li,wp,ws,wi,axp,axs,axi),intype)

#%% Kolonderski
def func_Kolonderski(ksx,ksy,kix,kiy,lp,ls,li,wp,ws,wi,axp,axs,axi):
    omega_p = lambda2omega(lp)
    omega_s = lambda2omega(ls)
    omega_i = lambda2omega(li)

    kpx = ksx + kix
    kpy = ksy + kiy
    kpz = np.sqrt(4*np.pi/lp**2 - kpx**2 - kpy**2)
    ksz = np.sqrt(4*np.pi/ls**2 - ksx**2 - ksy**2)
    kiz = np.sqrt(4*np.pi/li**2 - kix**2 - kiy**2)
    
    return wp*ws*wi*np.exp(-ws**2/2*(ksx**2 + ksy**2) - wi**2/2*(kix**2+kiy**2) - wp**2/2*(kpx**2+kpy**2)) * np.sinc(L/2*(kpz - ksz - kiz + 2*np.pi/Lambda))


def psi_Kolonderski(lp,ls,li,wp,ws,wi,axp,axs,axi):
# Integrale come somma; non funziona
#    (ksx,dksx) = (ksy,dksy) = (kix,dkix) = (kiy,dkiy) = np.linspace(-1000000,1000000,50,retstep=True)
#    
#    f = func_Kolonderski(ksx[:,None,None,None],ksy[None,:,None,None],kix[None,None,:,None],kiy[None,None,None,:],lp,ls,li,wp,ws,wi,axp,axs,axi)
#    return np.sum(np.sum(np.sum(np.sum(f))))*dksx*dksy*dkix*dkiy
# Integrale usando le funzioni di scipy: eterno
    return integ.nquad(func_Kolonderski,[[-np.inf,np.inf],[-np.inf,np.inf],[-np.inf,np.inf],[-np.inf,np.inf]],args=(lp,ls,li,wp,ws,wi,axp,axs,axi))

#%% Kolonderski paraxial
def func_Kolonderski(ksx,ksy,kix,kiy,lp,ls,li,wp,ws,wi,axp,axs,axi):
    omega_p = lambda2omega(lp)
    omega_s = lambda2omega(ls)
    omega_i = lambda2omega(li)

    kpx = ksx + kix
    kpy = ksy + kiy
    kpz = np.sqrt(4*np.pi/lp**2 - kpx**2 - kpy**2)
    ksz = np.sqrt(4*np.pi/ls**2 - ksx**2 - ksy**2)
    kiz = np.sqrt(4*np.pi/li**2 - kix**2 - kiy**2)
    
    return wp*ws*wi*np.exp(-ws**2/2*(ksx**2 + ksy**2) - wi**2/2*(kix**2+kiy**2) - wp**2/2*(kpx**2+kpy**2)) * np.sinc(L/2*(kpz - ksz - kiz + 2*np.pi/Lambda))


def psi_Kolonderski(lp,ls,li,wp,ws,wi,axp,axs,axi):
# Integrale come somma; non funziona
#    (ksx,dksx) = (ksy,dksy) = (kix,dkix) = (kiy,dkiy) = np.linspace(-1000000,1000000,50,retstep=True)
#    
#    f = func_Kolonderski(ksx[:,None,None,None],ksy[None,:,None,None],kix[None,None,:,None],kiy[None,None,None,:],lp,ls,li,wp,ws,wi,axp,axs,axi)
#    return np.sum(np.sum(np.sum(np.sum(f))))*dksx*dksy*dkix*dkiy
# Integrale usando le funzioni di scipy: eterno
    return integ.nquad(func_Kolonderski,[[-np.inf,np.inf],[-np.inf,np.inf],[-np.inf,np.inf],[-np.inf,np.inf]],args=(lp,ls,li,wp,ws,wi,axp,axs,axi))


#%% Phasematching
li = ls = np.arange(809.9e-9,810.1e-9,1e-12)
#psi1 = np.zeros((ls.size,li.size))
psi = np.zeros((ls.size,li.size))

#t0 = time.time()
#for i in range(ls.size):
#    for j in range(li.size):
#        psi[j,i] = np.abs(psi_Bennink(1/(1/ls[i]+1/li[j]),ls[i],li[j],wp,ws,wi,axp,axs,axi,'stupid'))
#t1 = time.time()

t2 = time.time()
def int_cycle(i):
    loc = np.zeros(li.size)
    for j in range(li.size):
        loc[j] = np.abs(psi_Bennink(1/(1/ls[i]+1/li[j]),ls[i],li[j],wp,ws,wi,axp,axs,axi,'stupid'))
    return loc

r = Parallel(n_jobs=8)(delayed(int_cycle)(i) for i in range(ls.size))

for i in range(ls.size):
    psi[:,i] = r[i]
t3 = time.time()

#print('Sequential in',t1-t0,'s')
print('Parallel in',t3-t2,'s')
   
#plt.figure()
#plt.plot(ls*1e9,psi1)
#plt.figure()
#plt.plot(ls*1e9,psi2)
#plt.figure()
#plt.plot(ls*1e9,psi3)
    
#%% Figures
plt.figure()
plt.imshow(psi, origin='lower', extent=(ls.min()*1e9, ls.max()*1e9, li.min()*1e9, li.max()*1e9),interpolation='nearest', cmap=cm.gist_rainbow)
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


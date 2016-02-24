# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 14:55:43 2016

@author: schiavon
"""

from util import lambda2omega
import numpy as np


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

#%% Kolonderski
def func_Kolonderski_full(ksx,ksy,kix,kiy,lp,ls,li,wp,ws,wi,axp,axs,axi):
    omega_p = lambda2omega(lp)
    omega_s = lambda2omega(ls)
    omega_i = lambda2omega(li)

    kpx = ksx + kix
    kpy = ksy + kiy
    kpz = np.sqrt(4*np.pi/lp**2 - kpx**2 - kpy**2)
    ksz = np.sqrt(4*np.pi/ls**2 - ksx**2 - ksy**2)
    kiz = np.sqrt(4*np.pi/li**2 - kix**2 - kiy**2)
    
    return wp*ws*wi*np.exp(-ws**2/2*(ksx**2 + ksy**2) - wi**2/2*(kix**2+kiy**2) - wp**2/2*(kpx**2+kpy**2)) * np.sinc(L/2*(kpz - ksz - kiz + 2*np.pi/Lambda))


def psi_Kolonderski_full(lp,ls,li,wp,ws,wi,axp,axs,axi):
# Integrale come somma; non funziona
#    (ksx,dksx) = (ksy,dksy) = (kix,dkix) = (kiy,dkiy) = np.linspace(-1000000,1000000,50,retstep=True)
#    
#    f = func_Kolonderski(ksx[:,None,None,None],ksy[None,:,None,None],kix[None,None,:,None],kiy[None,None,None,:],lp,ls,li,wp,ws,wi,axp,axs,axi)
#    return np.sum(np.sum(np.sum(np.sum(f))))*dksx*dksy*dkix*dkiy
# Integrale usando le funzioni di scipy: eterno
    return integ.nquad(func_Kolonderski,[[-np.inf,np.inf],[-np.inf,np.inf],[-np.inf,np.inf],[-np.inf,np.inf]],args=(lp,ls,li,wp,ws,wi,axp,axs,axi))

# Kolonderski paraxial
#def func_Kolonderski(lp,ls,li,wp,ws,wi,axp,axs,axi)

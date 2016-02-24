# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 14:55:23 2016

@author: schiavon
"""

from util import lambda2omega
import numpy as np
from myintegrate import complex_integral
from sellmeier import n
from scipy.constants import c

def func(z,lp,ls,li,wp,ws,wi,axp,axs,axi,Lambda):
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

def psi(lp,ls,li,wp,ws,wi,axp,axs,axi,L,Lambda,intype):
    return complex_integral(func,-L/2,L/2,(lp,ls,li,wp,ws,wi,axp,axs,axi,Lambda))
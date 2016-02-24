# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 14:54:48 2016

@author: schiavon
"""

import numpy as np
import scipy.integrate as integ

N = 10000

#def complex_integral(func,a,b,args,intype='stupid'):
#    real_func = lambda z: np.real(func(z,*args))
#    imag_func = lambda z: np.imag(func(z,*args))
#    if intype == 'quad':
#        real_int = integ.quad(real_func,a,b)
#        imag_int = integ.quad(imag_func,a,b)
##        print(real_int)
##        print(imag_int)
#        return real_int[0] + 1j * imag_int[0]
#    elif intype == 'quadrature':
#        real_int = integ.quadrature(real_func,a,b)
#        imag_int = integ.quadrature(imag_func,a,b)
##        print(real_int)
##        print(imag_int)
#        return real_int[0] + 1j * imag_int[0]
#    elif intype == 'romberg':
#        real_int = integ.romberg(real_func,a,b)
#        imag_int = integ.romberg(imag_func,a,b)
##        print(real_int)
##        print(imag_int)
#        return real_int + 1j * imag_int
#    elif intype == 'stupid':
#        Npoints = 10000
#        z,dz = np.linspace(a,b,Npoints,retstep=True)
#        real_int = np.sum(real_func(z))*dz
#        imag_int = np.sum(imag_func(z))*dz
##        real_int = integ.romb(real_func(z),dz)
##        imag_int = integ.romb(imag_func(z),dz)
##        print(real_int)
##        print(imag_int)
#        return real_int + 1j*imag_int

def complex_integral(func,a,b,args):
    real_func = lambda z: np.real(func(z,*args))
    imag_func = lambda z: np.imag(func(z,*args))
    Npoints = N
    z,dz = np.linspace(a,b,Npoints,retstep=True)
    real_int = np.sum(real_func(z))*dz
    imag_int = np.sum(imag_func(z))*dz
#        real_int = integ.romb(real_func(z),dz)
#        imag_int = integ.romb(imag_func(z),dz)
#        print(real_int)
#        print(imag_int)
    return real_int + 1j*imag_int

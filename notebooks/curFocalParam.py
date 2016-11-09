# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 16:23:14 2016

@author: schiavon
"""

# script to calculate the current focal paramters of the SPDC source

import numpy as np
import sellmeier as ppktp

def xi(L,wavelength,n,W):
    return L*wavelength/(2*np.pi*n*W**2)
    
    
# PPKTP crystal paramters
L = 0.03
Lambda = 10e-6

# beam from the CFC-11X-A collimator
W_CFC_0 = 770e-6
lp = 404.5e-9
#z0_CFC_0 = np.pi*W_CFC_0**2/lp

# lens LA1484-A
#n_NBK7 = 1.5303
#R1 = 0.1545 # m
#f1 = R1/(n_NBK7 - 1)
#f1 = 0.333

#print('Pump lens (LA1484-A), f =',f1)
#print('z0 from CFC-11X-A =',z0_CFC_0)

# waist at the crystal
#M1 = f1/z0_CFC_0
#Wp = W_CFC_0*M1
Wp = 53e-6

# parameters of the pump field
polaxis = 'y'
T = 20
n_p = ppktp.n(lp,polaxis,T)
xip = xi(L,lp,n_p,Wp)

print('Pump waist at the crystal =',Wp)
print('Xi_pump =',xip)


# idler and signal parameters
ls = li = 809e-9

# fiber parameters (from rp-photonics.com/mode_radius.html)
a = 2.2e-6 # core of the fibre
NA = 0.13 # numerical aperture
V = 2*np.pi*a*NA/ls
MFD = 2*a*(0.65 + 1.619/V**(3/2) + 2.879/V**6 - (0.016 + 1.561/V**7))
#MFD = 5e-6
W0s_fiber = MFD/2

# aspheric lens used (C280TME)
fasph = 18.4e-3

# waist at the crystal (signal and idler)
z_i = z_s = 0.245 # distance aspheric lens-crystal
z0_i = z0_s = np.pi*W0s_fiber**2/ls

Mi = Ms = np.sqrt( (z_s-fasph)**2/(fasph**2-z0_s**2) )

Wi = Ws = W0s_fiber*Ms

axis_signal = 'y'
axis_idler = 'z'

n_s = ppktp.n(ls,axis_signal,T)
n_i = ppktp.n(li,axis_idler,T)

print('Signal (idler) waist at the crystal =',Ws)
print('Xi_signal =',xi(L,ls,n_s,Ws))
print('Xi_idler =',xi(L,li,n_i,Ws))
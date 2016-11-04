# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 15:06:17 2016

@author: schiavon
"""

""" This script calculates the optimal lens distribution for Alice's and Bob's
    channels (in the hypothesis of collimation and focusing). """
    
import numpy as np
#import matplotlib.pyplot as plt
from scipy.constants import pi
import sellmeier as ppktp

from gaussbeam import Lens

n = 1.5106

f_asph = (8.011e-3,18.4e-3)
model = ('C240TME','C280TME')
wavelength = 809e-9

# fiber parameters (from rp-photonics.com/mode_radius.html)
a = 2.2e-6 # core of the fibre
NA = 0.13 # numerical aperture
V = 2*pi*a*NA/wavelength
MFD = 2*a*(0.65 + 1.619/V**(3/2) + 2.879/V**6 - (0.016 + 1.561/V**7))
#MFD = 5e-6
W0 = MFD/2

# due lenti e collimazione
print('Due lenti +  collimazione')
R = np.array([25.8e-3,38.6e-3,51.5e-3,64.4e-3,103e-3,128.8e-3,154.5e-3,257.5e-3,386.3e-3,515.1e-3])
f = ((n-1)/R)**-1
nominal_f = (50,75,100,125,200,250,300,500,750,1000)

for i in range(len(f_asph)):
    lns = Lens(f_asph[i],wavelength)
    z = f_asph[i]
    
    W01,z1 = lns.propagate(W0,z)
    print('Collimator:',model[i])
    print('New waist',W01*1e6,'um at z =',z1,'mm, with z0 =',W01**2*pi/wavelength*1e3,'mm')
    print()
    
    for j in range(len(f)):
        lns1 = Lens(f[j],wavelength)
        zsx = 0.3
        W02,z2 = lns1.propagate(W01,zsx)
        print('* lens of focal length',nominal_f[j],'mm')
        print('Waist at crystal',W02*1e6,'um at z =',z2,'mm from the lens')
        print('Bennink xi =',0.03*wavelength/(2*pi*ppktp.n(wavelength,'z',T=20)*W02**2))
        print()
    print('\n')
    
#%% singola lente
print('Single aspheric lens')
print('MFD =',MFD*1e6,'um')

W01 = 47e-6
W01 = 30e-6
print('Bennink xi =',0.03*wavelength/(2*pi*ppktp.n(wavelength,'y',T=20)*W01**2))

M = W01/W0
z0 = pi*W0**2/wavelength

for i in range(len(f_asph)):
    z = f_asph[i] + np.sqrt( (f_asph[i]/M)**2 - z0**2 )
    z1 = f_asph[i] + M**2*(z-f_asph[i])
    print('Collimator:',model[i])
    print('Distance lens-fiber:',z*1e3,'mm')
    print('Distance lens-crystal:',z1*1e3,'mm')
    print('Beam radius at lens:',W0*np.sqrt(1+(z/z0)**2)*1e3,'mm with lens radius',5.50/2,'mm')
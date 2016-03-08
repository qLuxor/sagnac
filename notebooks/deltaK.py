# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 17:33:38 2016

@author: schiavon
"""

""" This script shows the variation of the phase mismatch as a function of the 
    temperature for a PPKTP with 10 um grating period."""
    
import numpy as np
import matplotlib.pyplot as plt
from sellmeier import n
from scipy.constants import pi

lp = 405e-9
ls = li = 810e-9

# take into account the thermal dilation of the crystal
Lambda0 = 10e-6
Lambda = lambda T: Lambda0*(1 + 7e-6*(T-25) + 4.4e-9*(T-25)**2)

axp = axi = 'y'
axs = 'z'

T = np.arange(10,50,0.0001)

Dk = 2*pi*(n(lp,axp,T)/lp - n(ls,axs,T)/ls - n(li,axi,T)/li - 1/Lambda(T))

plt.plot(T,Dk)
plt.ylabel(r'Phase mismatch $\Delta k$')
plt.xlabel(r'Temperature [${}^\circ$ C]')
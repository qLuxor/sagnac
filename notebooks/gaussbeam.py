# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 16:12:25 2016

@author: schiavon
"""

from scipy.constants import pi
import numpy as np

class Lens:
    
    def __init__(self,focus,wavelength):
        self.f = focus
        self.l = wavelength
        
    def propagate(self,W0,z):
        z0 = pi*W0**2/self.l
        M = self.f/np.sqrt((z-self.f)**2 + z0**2)
        W01 = W0*M
        z1 = self.f + M**2*(z-self.f)
        return W01,z1
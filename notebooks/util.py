# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 15:02:17 2016

@author: schiavon
"""

import numpy as np
from scipy.constants import c


def omega2lambda(omega):
    return 2*np.pi*c/omega

def lambda2omega(l):
    return 2*np.pi*c/l

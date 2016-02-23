# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 12:49:31 2015

@author: schiavon
"""

import numpy as np

def sellmeier_ktp(wl,axis):
    if axis=='x':
        res = np.sqrt(3.0065 + 0.03901/((wl*1e6)**2 - 0.04251) - 0.01327*(wl*1e6)**2)
    if axis=='y':
        res = np.sqrt(3.0333 + 0.04154/((wl*1e6)**2 - 0.04547) - 0.11408*(wl*1e6)**2)
    if axis=='z':
        res = np.sqrt(3.3134 + 0.05694/((wl*1e6)**2 - 0.05658) - 0.01682*(wl*1e6)**2)
    return res

def t_correction_ktp(wl,T,axis):
    if axis=='y':
        n_1 = 6.2897e-6 + 6.3061e-6/(wl*1e6) - 6.0629e-6/((wl*1e6)**2) + 2.6486e-6/((wl*1e6)**3)
        n_2 = -0.14445e-8 + 2.2244e-8/(wl*1e6) - 3.5770e-8/((wl*1e6)**2) + 1.3470/((wl*1e6)**3)
        res = n_1 * (T-20) + n_2 * (T-20)**2
    elif axis=='z':
        n_1 = 9.9587e-6 + 9.9228e-6/(wl*1e6) - 8.9603e-6/((wl*1e6)**2) + 4.1010e-6/((wl*1e6)**3)
        n_2 = -1.1882e-8 + 10.459e-8/(wl*1e6) - 9.8138e-8/((wl*1e6)**2) + 3.1481/((wl*1e6)**3)
        res = n_1 * (T-20) + n_2 * (T-20)**2
    return res
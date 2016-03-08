# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 12:49:31 2015

@author: schiavon
"""

import numpy as np

def sellmeier_ktp(wl,axis):
    ''' Sellmeier equations for PPKTP ('waveguide nonlinear optic devices', at 20°C),
        x and y are ordinary axes, while z is the extraordinary one.
        l is the wavelength in m'''
    if axis=='x':
        res = np.sqrt(3.0065 + 0.03901/((wl*1e6)**2 - 0.04251) - 0.01327*(wl*1e6)**2)
    if axis=='y':
        res = np.sqrt(3.0333 + 0.04154/((wl*1e6)**2 - 0.04547) - 0.11408*(wl*1e6)**2)
    if axis=='z':
        res = np.sqrt(3.3134 + 0.05694/((wl*1e6)**2 - 0.05658) - 0.01682*(wl*1e6)**2)
    return res

def t_correction_ktp(wl,T,axis):
    ''' Temperature correction for Sellmeier equations (from Emanueli et. al.,
        Applied Optics Vol. 42, No 33).
        The function return 0 if axis=='x'.'''
    if axis=='y':
        n_1 = 6.2897e-6 + 6.3061e-6/(wl*1e6) - 6.0629e-6/((wl*1e6)**2) + 2.6486e-6/((wl*1e6)**3)
        n_2 = -0.14445e-8 + 2.2244e-8/(wl*1e6) - 3.5770e-8/((wl*1e6)**2) + 1.3470/((wl*1e6)**3)
        res = n_1 * (T-20) + n_2 * (T-20)**2
    elif axis=='z':
        n_1 = 9.9587e-6 + 9.9228e-6/(wl*1e6) - 8.9603e-6/((wl*1e6)**2) + 4.1010e-6/((wl*1e6)**3)
        n_2 = -1.1882e-8 + 10.459e-8/(wl*1e6) - 9.8138e-8/((wl*1e6)**2) + 3.1481/((wl*1e6)**3)
        res = n_1 * (T-20) + n_2 * (T-20)**2
    else:
        res = 0
    return res


def n_Kato(l,axis):
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


def n(l,axis,T=20,typ='christ'):
    ''' Returns the value of the refractive index along the axis at wavelength l
        and temperature T, using the equations derived from Kato-Takaoka (if 
        typ == 'kato') or those found in the thesis of A. Christ (if typ == 'christ').'''
    if typ == 'christ':
        return sellmeier_ktp(l,axis) + t_correction_ktp(l,T,axis)
    else:
        return n_Kato(l,axis) + t_correction_ktp(l,T,axis)
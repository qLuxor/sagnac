# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 18:37:10 2016

@author: schiavon
"""

# herlding efficiency of SPDC (directly to detectors)

import numpy as np

hdA = hdB = 0.6
hfA = 1.77/1.85
hfB = 1.81/1.87

ch1 = np.array([775,865,332,366])
ch2 = np.array([792,895,337,380])
c12 = np.array([109,139,45,58])
c12_acc = np.array([1.6,1.5,0.3,0.33])

eta = (c12-c12_acc)/np.sqrt(ch1*ch2)

hA = np.mean((c12-c12_acc)/ch1/hdA)
hB = np.mean((c12-c12_acc)/ch1/hdB)

print('eta =',np.mean(eta))
print('hA =',np.mean(hA))
print('hB =',np.mean(hB))

c12_40 = np.mean(c12[2:3])
ch1_40 = np.mean(ch1[2:3])
ch2_40 = np.mean(ch2[2:3])

c12_60 = np.mean(c12[0:1])
ch1_60 = np.mean(ch1[0:1])
ch2_60 = np.mean(ch2[0:1])

I_40 = np.mean([2.1,2.6])
I_60 = np.mean([5.7,6.2])

print('Coinc_40 =',c12_40/I_40,'kcoinc/mW')
print('Coinc_60 =',c12_60/I_60,'kcoinc/mW')

dl = 0.18

print('Production rate I40 =',c12_40/hdA/hdB/I_40,'kHz/mW')
print('Production rate I60 =',c12_60/hdA/hdB/I_60,'kHz/mW')

print('Production rate I40 =',c12_40/hdA/hdB/I_40/dl,'kHz/mW/nm')
print('Production rate I60 =',c12_60/hdA/hdB/I_60/dl,'kHz/mW/nm')
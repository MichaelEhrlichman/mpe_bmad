#!/bin/env python3

import numpy as np

c = 299792458.
me = 510999.

laser = 1030e-9
E = 98e6

g = 0.032  #undulator gap in meters
lu = 0.054  #undulator period

Bu = 2.2 * np.exp(-np.pi*g/lu) * np.sin(np.pi/4.0) / (np.pi/4.0)

K = Bu * lu * c / 2 / np.pi / me

print("Bu = {} T".format(Bu))
print("K = {}".format(K))


gami = E/me
K_formula = 2*np.sqrt(laser*gami**2/lu-0.5)
B_formula = 2.0*np.pi*me/c * K_formula / lu
g_formula = -lu/np.pi*np.log(np.pi/4.0*B_formula/2.2/np.sin(np.pi/4.0))

print("K from formula = {}".format(K_formula))
print("B from formula = {}".format(B_formula))
print("g from formula = {}".format(g_formula))

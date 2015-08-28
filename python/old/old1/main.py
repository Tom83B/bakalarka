# -*- coding: utf-8 -*-
import os.path
import sys
sys.path.append('/home/tom/Documents/skola/bakalarka/python/')
sys.path.append('/home/tom/Documents/skola/bakalarka/python/EOS/')
import pdb

import matplotlib.pyplot as plt
import numpy as np

#nahrání veličin pro okolí a opacitní tabulky
import ext_data as ext
from ext_data import g,r,z
from const import *

#jen do hloubky cca 12Mm
z = ext.z[0:679]

#výpočet pozaďového statického řešení
execfile('/home/tom/Documents/skola/bakalarka/python/static2.py')

def pder(p, z):
	r = rho(Tfunc(z), np.log(p), z)
	D = r*g(z) 
	return D

pi0 = pe[499]-B0*B0/(8e-7*np.pi)
zinv = z[::-1]
pi = integrate.odeint(pder, pi0, zinv)
pi = pi[::-1]

y0 = np.linspace(np.sqrt(Btop), np.sqrt(B0), 500)

#výpočet magnetického pole
execfile('/home/tom/Documents/skola/bakalarka/python/magnetic.py')

Y = comp_mag(y0)

#funkce pro uložení fcí obsahujících opacitu
execfile('/home/tom/Documents/skola/bakalarka/python/plots.py')

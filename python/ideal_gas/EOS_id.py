#coding: utf-8 -*-
import os.path
import sys
sys.path.append('/home/tom/Documents/skola/bakalarka/python/')
sys.path.append('/home/tom/Documents/skola/bakalarka/python/EOS/')

import numpy as np
from scipy.optimize import fixed_point
from const import *

import ext_data as ext

#ionizace
def x(T, P):
	A = 0.5*(2./P)*((np.sqrt(2*np.pi*me)/planck)**3)*((k*T)**(5./2.))*np.exp(-chi/(k*T))
	return np.sqrt(A/(1+A))

#molární hmotnost opravená o ionizaci
def mu(T, P):
	return mu0/(1+x(T, P))

def opac(T, P, z):
    opac0 = ext.opac(T, rho(T,P,z), z)
    return opac0[0,0]

def rho(T, P, z):
	return 1e-3*P*mu(T,P)/(R*T)

def PofRho(T, rho):
	return fixed_point(lambda P: 1e3*rho*R*T/mu(T,P), 1e3*rho*R*T)

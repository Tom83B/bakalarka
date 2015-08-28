# -*- coding: utf-8 -*-
import numpy as np
from scipy import integrate
from scipy.optimize import bisect

from common.const import Btop, alpha
from common.derivatives import derivative 
import loadbcg as bcg
from magnetic import comp_mag
from EOS import adG, S, TofS, opac, rho, PofRho
import EOS
from ext_data import g
import ext_data as ext

import pdb

#    atmosfera
#____________________________________________________________

top = 0
T_spot = 3500.

T_eff = T_spot*2**0.25
S0 = S(bcg.T[bcg.z.size-1], bcg.P[bcg.z.size-1], 12.5e6)
print S0
tau_top = 2./3.
tau = np.linspace(0, tau_top, 50)
T_atm = T_eff*(0.5*(1+tau*3./2.))**0.25
T_atm = interp(tau, T_atm)
print T_atm(tau_top)

def atm_der(y,tau):
	z, P = y
	return np.array([ 1./(EOS.opac(T_atm(tau), P, z)*EOS.rho(T_atm(tau), P, z)), g(z)/EOS.opac(T_atm(tau), P, z) ])

def atm_entrpy(tau0):
	tau = np.linspace(tau0, tau_top, 50)
	yinit = [ top, P0 ]

	atm = integrate.odeint(atm_der, yinit, tau)
	z_atm = atm[:,0]
	P_atm = atm[:,1]
	return S(T_atm(tau_top), max(P_atm), max(z_atm))
	
	yinit = np.array([ T_atm(tau_top), np.log(P_atm[P_atm.size-1]) ])
	z = np.linspace(z0, 12.5e6, 500)
	y = integrate.odeint(derivative, yinit, z)
	T = y[:,0]
	P = np.exp(y[:,1])
	print S(max(T), max(P), max(z))
	return S(max(T), max(P), max(z))

P0 = bcg.P[0]*(1-(bcg.T[0]-T_atm(0))/bcg.T[0])
#tau0 = bisect(lambda P: atmosphere_enthropy(P)-S0, 10000., 30000.)
P0 = EOS.PofRho(T_spot,1e-5)

yinit = np.array([ top, P0 ])
atm = integrate.odeint(atm_der, yinit, tau)
z_atm = atm[:,0]
P_atm = atm[:,1]
print z_atm

top = max(z_atm)
top = 0

#____________________________________________________________

z = np.linspace(top, 12.5e6, 500)

def surface_temperature(B):
	dP = B*B/(8e-7*np.pi)
	P0 = bcg.P[499]-dP
	T0 = TofS(P0, S(bcg.T[499],bcg.P[499],z[499]),bcg.T[499],z[499])
	yinit = np.array([ T0, np.log(P0) ])
	z_inv = np.linspace(12.5e6, top, 500)
	a = integrate.odeint(derivative, yinit, z_inv)
	print B, a[499,0]
	return a[499,0]

#B0 = bisect(lambda B: surface_temperature(B)-T_atm(tau_top), 9., 11.) #treba nastavit vysokou toleranci ~0.1K
B0 = 9.83040618896
dP = B0*B0/(8e-7*np.pi)
#P0 = bcg.P[499]-dP
#T0 = TofS(P0, S(bcg.T[499],bcg.P[499],z[499]),bcg.T[499],z[499])
T0 = T_atm(tau_top)
T0 = 3500.
P0 = max(P_atm)
#P0 = bisect(lambda P: S(T0,P,0)-S0,1e1,bcg.P[0])
P0 = bcg.P[0]*(T0/bcg.T[0])**2.5

yinit = np.array([ T0, np.log(P0) ])
z_inv = np.linspace(12.5e6, top, 500)
a = integrate.odeint(derivative, yinit, z)
T = a[:,0]
#T = T[::-1]
P = np.exp(a[:,1])
#P = P[::-1]
#z = z_inv[::-1]
print S(max(T), max(P), max(z))

#z = np.concatenate([ z_atm, z[1:]] )
#P = np.concatenate([ P_atm, P[1:] ])
#T = np.concatenate([ T_atm(tau), T[1:] ])

#pdb.set_trace()

#y0 = np.linspace(np.sqrt(Btop), np.sqrt(B0), 500)

#y = comp_mag(y0, P)

#pdb.set_trace()

data = np.concatenate((np.array([z]).T, np.array([T]).T, np.array([np.exp(lnP)]).T), axis=1)
np.savetxt('initial_state.dat', data)

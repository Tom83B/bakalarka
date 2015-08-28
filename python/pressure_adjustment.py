import numpy as np
from scipy import optimize, integrate

from const import phi
from EOS import rho

def mass(P0, y, z, T):
	P = integrate.odeint(PDeriv, P0, np.linspace(1.25e7, 0, 500), args = (T))
	y = comp_mag(y, P)
	B = y*y

	rhoArr = np.zeros(500)
	for i in range(500):
		rhoArr[i] = rho(T[i], P[i], z[i])
	mmbrs = np.zeros(499)
	for i in range(499):
		mmbrs[i] = (rhoArr[i]/B[i]+rhoArr[i+1]/B[i+1])*(z[i+1]-z[i])/2
	integral = mmbrs.sum()
	return phi*integral

def adjust_baseP(y, z, T, rho0, known_mass = 0, static={"known_mass":0}):
	return optimize.newton(lambda P: mass(P, y, z, T)-static["known_mass"]-phi*rho0*v0*dt/(y[499]*y[499]))

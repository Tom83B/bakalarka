# -*- coding: utf-8 -*-
import os.path
import sys
sys.path.append('/home/tom/Documents/skola/bakalarka/python/')

from derivatives import Frad, Fconv
from scipy.interpolate import UnivariateSpline
from EOS import rho, cp
from const import *
from ext_data import opac


for i in range (0, z.size):
	FconvArr[i] = Fconv(T[i], P[i], z[i])
	FradArr[i] = Frad(T[i], P[i], z[i])
	rhoArr[i] = rho(T[i], P[i], z[i])
	cpArr[i] = cp(T[i], P[i])
	opacArr[i] = opac(T[i],rhoArr[i])

F = FconvArr+FradArr
F_spline = UnivariateSpline(z, F, s=0)
dFdz = F_spline.derivative()
dFdz = dFdz(z)

#dFdz = np.diff(F)/h
#dFdz = np.append(dFdz, dFdz[dFdz.size-1])

#  semi-implicitni metoda pro casovy prirustek

#d_main = np.empty(500)
#d_sub = np.zeros(500)
#d_sup = np.zeros(500)
#b = np.zeros(500)
#
#h = 25000	#grid spacing
#
#Bz0 = np.sqrt(8e-7*np.pi*(bcg.P[499]-P[499]))
#b[0] = np.sqrt(Btop)
#b[499] = np.sqrt(Bz0)
#
#for i in range(1,499): d_main[i] = rhoArr[i]*cpArr[i]+(dt/h)
#for i in range(0,498): d_sub[i] = -(eps*phi*y[i+1])/(2*np.pi*h*h)
#for i in range(2,500): d_sup[i] = -(eps*phi*y[i-1])/(2*np.pi*h*h)
#for i in range(1,499): b[i] = y[i]**4+(eps*8e-7*np.pi*(bcg.P[i]-P[i]))
#d_main[0] = 1
#d_main[499] = 1
#
#data = [d_sub, d_main, d_sup]
#diags = [-1,0,1]
#A = sparse.spdiags(data, diags, 500, 500, format='csc')
#
#return spsolve(A,b)
#
#dTdt = -dFdz/(cpArr*rhoArr)
#T = T+dt*dTdt
#P0 = adjust_baseP(y, z, T, rho(T[499], P[499], z[499]))
#P = integrate.odeint(lambda P, z: PDeriv(T_func(z), P, z), P0, z_inv)

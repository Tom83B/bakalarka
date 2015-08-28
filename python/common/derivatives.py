# -*- coding: utf-8 -*-
import os.path
import sys
sys.path.append('/home/tom/Documents/skola/bakalarka/python/')

import numpy as np
from scipy.optimize import newton
from const import *
from EOS import U, V, radG, adG, rho, cp, H, mu, delta
from OPAL.opacity import opac
from ext_data import g,M
import ext_data as ext
#from xi_interp import xi

#   derivace teploty
#__________________________________________________________________

#def messyTDeriv(U, V, radG, adG):
#	TLDR = V*(1./(3*V)-(2**(1./3.)*(-1-6.*V))/(3.*V*(2.+18.*V+27.*adG*V**2+27*radG*V**2+np.sqrt(4*(-1-6*V)**3+(2+18*V+27*adG*V**2+27*radG*V**2)**2))**(1./3.))+(1./(3.*2**(1./3)*V))*((2+18*V+27*adG*V**2+27*radG*V**2+np.sqrt(4*(-1-6*V)**3+(2+18*V+27*adG*V**2+27*radG*V**2)**2))**(1./3.)))**3.+radG
#	return TLDR

def messyTDeriv(U, radG, adG):
	W = radG-adG
    #A = (19*U)/27-(184*2**(1./3.)*U**2)/(27*(1168*U**3+2187*U*W+27*np.sqrt(3)*np.sqrt(2048*U**6+2336*U**4*W+2187*U**2*W**2))**(1./3.))+1/27*2**(2./3.)*(1168*U**3+2187*U*W+27*np.sqrt(3)*np.sqrt(2048*U**6+2336*U**4*W+2187*U**2*W**2))**(1./3.)
    #A = xi(U,W)
	A = newton(lambda xi: (xi-U)**3+(8*U/9)*(xi**2-U**2-W), adG, lambda xi: 3*(xi-U)**2+2*8*U*xi/9)
#	A = newton(lambda xi: (xi-U)**3+(f/b)*U*(xi**2-U**2-W), (adG+radG)/2, lambda xi: 3*(xi-U)**2+2*(f/b)*U*xi)
    #print U, W, A
	return A**2-U**2+adG

def TDeriv(T, P, z):
	rG = radG(T, P, z)
	aG = adG(T, P)
#	aG = ext.adG_T(T)
#	cG = messyTDeriv( U(T,P,z), V(T,P,z), rG, aG )
	#print z, rG, aG, dtdz
	if rG<aG:
		return rG
	else:
		cG = messyTDeriv( U(T,P,z), rG, aG )
		return cG
#		return aG

#def TDeriv(T, P, z):
#	return max(radG(T, P, z), adG(T, P))

#   derivace tlaku
#__________________________________________________________________

def PDeriv(P, z, T):
	r = rho(T, P, z)
	D = r*g(z) 
	return D

#    derivace pro integraci
#__________________________________________________________________

def logPderivative(y, lnP,):
	z, lnT = y
	T = np.exp(lnT)
	P = np.exp(lnP)
	return np.array([ H(T, P, z), TDeriv(T, P, z) ])


def derivative(y, z):
	T, lnP = y
	P = np.exp(lnP)
#	print z, T, P
#	return np.array([ T*TDeriv(T, P, z)/H(T, P, z), 1/H(T, P, z) ])
	return np.array([ T*TDeriv(T, P, z)/H(T, P, z), rho(T,P,z)*g(z)/P ])
#	return np.array([ T*adG(T, P)/H(T, P, z), rho(T,P,z)*g(z)/P ])
#	return np.array([ T*ext.adG(z)*H(T, P, z), rho(T,P,z)*g(z)/P ])
#	return np.array([ T*TDeriv(T, P, z)/H(T, P, z), (np.log(ext.p(z+1000))-np.log(ext.p(z-1000)))/2000. ])
#	return np.array([ T*ext.adG(z)/H(T, P, z), (np.log(ext.p(z+1000))-np.log(ext.p(z-1000)))/2000. ])
#	return np.array([ (ext.T(z+1000)-ext.T(z-1000))/2000., rho(T,P,z)*g(z)/P ])
#	return np.array([ (ext.T(z+100)-ext.T(z-100))/200., (np.log(ext.p(z+100))-np.log(ext.p(z-100)))/200. ])


#	toky energie
#_______________________________________________________

def flux(T, dTdz, P, dPdz, z, alpha):
	G = (P/T)*(dTdz/dPdz)
	rG = radG(T, P, z)
#	rho = dPdz/g(z)
	r = rho(T, P, z)
	o = opac(T, r)
	u = U(T,P,z)
	aG = adG(T,P)
#	frad = 16*sigma*T**3*g(z)/(3*opac(T, P, z)*r*P)*dTdz/dPdz
	frad = (16./3.)*sigma*g(z)*T**4*G/(o*P)
#	if frad>sigma*Teff**4:
#		return np.array([ sigma*Teff**4, 0, sigma*Teff**4 ])
#	else:
	eG = aG-2*u**2+2*u*np.sqrt(G-aG+u**2)
	fconv = r*cp(T, P, z)*T*np.sqrt(g(z)*delta(T,P,z))*alpha**2*np.sqrt(H(T,P,z))*(G-eG)**1.5/(4*np.sqrt(2.))
#	fconv = r*cp(T, P, z)*T*np.sqrt(g(z)*delta(T,P,z))*alpha**2*np.sqrt(H(T,P,z))*(abs(G-eG))**(3./2.)/(4*np.sqrt(2.))
#	return np.array([ frad, fconv ])
#	ftot = c*3.846e26*g(z)/(4*M(z)*G*np.pi)
	ftot = (16./3.)*sigma*g(z)*T**4*rG/(o*P)
	return np.array([ frad, fconv, ftot ])

def Frad(T, P, z):
	dTdz = T*TDeriv(T, P, z)/H(T, P, z)
	FR = -(16*sigma*T**3/(3*opac(T, P, z)*rho(T, P, z)))*dTdz
	return FR

def Fconv(T, P, z, alpha):
	u = U(T, P, z)
	logG = TDeriv(T, P, z)
	logadG = adG(T, P) #definovano v temp_derivative
	#commaG = logadG-2*u**2+2*u*np.sqrt(logG-logadG+u**2)
	return rho(T, P, z)*cp(T, P, z)*T*np.sqrt(g(z)*delta(T, P, z))*alpha**2*np.sqrt(H(T,P,z))*8*U(T,P,z)*(radG(T,P,z)-logG)/(4*9*np.sqrt(2.))
#	return -b*np.sqrt(a*R/mu(T, P))*np.sqrt(alpha)*rho(T, P, z)*T**(1.5)*(f/b)*(radG(T,P,z)-logG)

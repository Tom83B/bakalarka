#coding: utf-8 -*-
import os.path
import sys
sys.path.append('/home/tom/Documents/skola/bakalarka/python/')
sys.path.append('/home/tom/Documents/skola/bakalarka/python/EOS/')

import numpy as np
from scipy.optimize import newton, fixed_point, bisect

import eos
from const import *
import ext_data as ext
from ext_data import g,M
import EOS_id

irad = 0

def S(T, P, z):
	reteos = eos.esac(X, ztab, T*1e-6, rho(T, P, 0)*1e-3, 9, irad)
	P, E, S, cv, chiRho, chiT, gamma1, nabla_ad_inv = reteos[:8]
	print P, S
	return S

def TofS(P, S0, T0, z):
	return newton(lambda T: S0-S(T,P,z), T0)

#ionizace
def x(T, P):
	A = 0.5*(2./P)*((np.sqrt(2*np.pi*me)/planck)**3)*((k*T)**(5./2.))*np.exp(-chi/(k*T))
	return np.sqrt(A/(1+A))

#molární hmotnost opravená o ionizaci
def mu(T, P):
	return mu0/(1+x(T, P))

def rho(T, P, z):
#	try:
	return eos.rhoofp(X, ztab, T*1e-6, P*1e-11, irad)*1e3*1
#	except:
#	return bisect(lambda r: EOS.PofRho(3500,r)-P, 1e-7, 10)
#	return eos.rhoofp(X, ztab, ext.T(z)*1e-6, ext.p(z)*1e-11, irad)*1e3*1
#	return ext.r(z)/1.5
#	return 1e-3*ext.p(z)*mu(ext.T(z),ext.p(z))/(R*ext.T(z))
#	return 1e-3*P*mu(T,P)/(R*T)

def opac(T, P, z):
    opac0 = ext.opac(T, rho(T,P,z), z)
    return opac0


#tepelná kapacita
def cp(T, P, z):
#	reteos = eos.esac(X, ztab, (T)*1e-6, rho(T, P, z)*1e-3, 9, irad)
#	P, E, S, cv, chiRho, chiT, gamma1, nabla_ad_inv = reteos[:8]
#	print P
#	c= 1e2*cv+(P/(T*rho(T,P,z)))*chiT*chiT/chiRho
#	return c
#	return 1e3*R/(mu(T, P))*(2.5+0.5*x(T, P)*(1-x(T, P))*(2.5+chi/(k*T))**2)
#	r = eos.rhoofp(X, ztab, T*1e-6, np.exp(lnP)*1e-11, irad)	#rhoofp vraci hustotu v g/cm^3
	h = 100.
#	print 'a'
	reteos_plus = eos.esac(X, ztab, (T+h)*1e-6, rho(T, P, z)*1e-3, 7, irad)
#	print 'b'
	reteos_minus = eos.esac(X, ztab, (T-h)*1e-6, rho(T, P, z)*1e-3, 7, irad)
	P, E_plus, S_plus, cv, chiRho, chiT, gamma1 = reteos_plus[:7]
	P, E_minus, S_minus, cv, chiRho, chiT, gamma1 = reteos_minus[:7]
	return 1e2*T*(1./(2*h))*(S_plus-S_minus)
#	return 1e8*(E_plus-E_minus)/(2*h)+P*delta(T,P,z)/(rho(T,P,z)*T)
#	def Q(Ta):
#		return k*Ta*(3./2.)/(mu(Ta,P)*m_U)
#	ret = (1/(2*h))*(Q(T+h)-Q(T-h))
#	return ret

#tlaková škála
def H(T, P, z):
	return P/(rho(T, P, z)*g(z))


#   gradienty
#_______________________________________________________

def V(T, P, z):
	V0 = (3./32)*(cp(T,P)*rho(T,P,z)*opac(T,P,z)*H(T,P,z))/(sigma*T**3)*alpha**2*np.sqrt(R*T/(8e-3*mu(T,P)))
	return V0

#def U(T, P, z):
#	U0 = 12*(sigma*T**3/(cp(T,P)*rho(T,P,z)**2*opac(T,P,z)*H(T,P,z)))*(1/alpha)**2*np.sqrt(8e-3*mu(T,P)/(R*T))
#	return U0


#def U(T, P, z):
#	U0 = 1/(f*np.sqrt(a))*12*(sigma*T**3/(cp(T,P)*rho(T,P,z)**2*opac(T,P,z)*H(T,P,z)**2))*(1/alpha)**2*np.sqrt(8*H(T,P,z)/(1*g(z)))
#	return U0

def delta(T, P, z):
#	h = 100.
#	r = rho(T, P, z)
#	r_plus = rho(T+h, P, z)
#	r_minus = rho(T-h, P, z)
#	drdt = (1/(2.*h))*(r_plus-r_minus)
#	return -(T/r)*drdt
	reteos = eos.esac(X, ztab, T*1e-6, rho(T, P, 0)*1e-3, 9, irad)
	P, E, S, cv, chiRho, chiT, gamma1, nabla_ad_inv = reteos[:8]
	return chiT/chiRho
#	return 1.

def PofRho(T,rho):
	reteos = eos.esac(X, ztab, T*1e-6, rho*1e-3, 9, irad)
	P, E, S, cv, chiRho, chiT, gamma1, nabla_ad_inv = reteos[:8]
	return P*1e11

def U(T, P, z):
	U0 = 12*(sigma*T**3/(cp(T,P,z)*rho(T,P,z)**2*opac(T,P,z)*H(T,P,z)**2))*(1/alpha)**2*np.sqrt(8*H(T,P,z)/(delta(T,P,z)*g(z)))
	return U0

def messyTDeriv(U, radG, adG):
	W = radG-adG
	A = newton(lambda xi: (xi-U)**3+(8*U/9)*(xi**2-U**2-W), (adG+radG)/2, lambda xi: 3*(xi-U)**2+2*8*U*xi/9)
	return A**2-U**2+adG

def TDeriv(T, P, z):
	return adG(T,P)

def radG(T, P, z):
#	Hi = (16./(3.*rho(T, P, z)*opac(T, P, z)))*sigma*T**4*(TDeriv(T, P, z)/H(T, P, z))
#	rG = 3*opac(T, P, z)*P*(4.*Hi)/(16*a*g(z)*c*T**4)
	rG = (3./64)*(opac(T,P,z)*P*3.846e26/(M(z)*np.pi*G*sigma*T**4))
	return rG

def adG(T, P):
	reteos = eos.esac(X, ztab, T*1e-6, rho(T, P, 0)*1e-3, 9, irad)
	P, E, S, cv, chiRho, chiT, gamma1, nabla_ad_inv = reteos[:8]
#	adG1 = rho(T, P, 0)
#	print 1/nabla_ad_inv-adG1
#	print adG0, adG1, rho(T,P,0), delta(T, P, 0)
	return 1/nabla_ad_inv
#	return adG1

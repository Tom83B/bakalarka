# -*- coding: utf-8 -*-
from scipy import integrate

#	definování veličin:
#_____________________________________________________

#tlaková škála
def H(T, lnP, z):
    return 1e3*(R*T)/(mu(T,lnP)*g(z))

#ionizace
def x(T, lnP):
	A = (2.5*np.log10(T)-13.53*5040/T-0.48-np.log(X)-lnP/np.log(10)-np.log10(mu0))
	return np.sqrt(10**A/(1+10**A))

#molární hmotnost opravená o ionizaci
def mu(T, lnP):
    return mu0/(1+x(T, lnP))

def rho(T, lnP,z):
    return 1e-3*mu(T,lnP)*np.exp(lnP)/(R*T)

def opac(T, lnP,z):
    opac0 = ext.opac(T, rho(T,lnP,z))
    return opac0[0,0]

#tepelná kapacita
def cp(T, lnP):
    return 1e3*R/(mu(T, lnP))*(2.5+0.5*x(T, lnP)*(1-x(T, lnP))*(2.5+chi/(k*T))**2)

#	toky energie
#_______________________________________________________

#načtení (d log T/d log p) + funkce navíc
execfile('temp_derivative.py')

def Frad(T,lnP,z):
	dTdz = T*TDeriv(T,lnP,z)/H(T, lnP, z)
	return  -(16*sigma*T**3/(3*opac(T,lnP,z)*rho(T,lnP,z)))*dTdz

def Fconv(T,lnP,z):
	u = 1/(f*np.sqrt(a))*alpha**(-2)*12*sigma*T**3/(cp(T,lnP)*rho(T,lnP,z)*opac(T,lnP,z)*H(T,lnP,z)**2)*np.sqrt(H(T,lnP,z)/g(z))
	logG = TDeriv(T, lnP, z)
	logadG = adG(T,lnP) #definovano v temp_derivative
	commaG = logadG-2*u**2+2*u*np.sqrt(logG-logadG+u**2)

	return -b*np.sqrt(a*R/mu(T,lnP))*np.sqrt(alpha)*rho(T,lnP,z)*T**(1.5)*(logG-commaG)**(1.5)

#	integrace rovnic
#____________________________________________________________

def derivative(y, z):
	T, lnP = y
#	return np.array([ T*TDeriv(T,lnP,z)/H(T, lnP, z), 1/H(T,lnP,z) ])
	return np.array([ T*TDeriv(T,lnP,z)/H(T, lnP, z), (np.log(ext.p(z+200))-np.log(ext.p(z)))/200. ])
#	return np.array([ (np.log(ext.p(z+25000))-np.log(ext.p(z)))/25000., T*TDeriv(T,lnP,z)/H(T, lnP, z) ])
#	return np.array([ (ext.T(z+100)-ext.T(z))/100., (np.log(ext.p(z+100))-np.log(ext.p(z)))/100. ])
#	return np.array([ (ext.T(z+100)-ext.T(z))/100., 1/H(T, lnP, z) ])


#počáteční podmínky
P0 = (2./3.)*g(0)/opac(ext.T(0), np.log(ext.p(0)),0)
P0 = ext.p(0)
Teff = 6000
Teff = ext.T(0)
yinit = np.array([ Teff, np.log(P0) ])

#integrace pomocí isoda pack
z = np.linspace(0., 2e7, 100)
y = integrate.odeint(derivative, yinit, z, mxstep=50000)
T = y[:,0]
lnP = y[:,1]

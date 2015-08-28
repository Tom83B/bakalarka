# -*- coding: utf-8 -*-
from scipy import integrate
from scipy import interpolate

irad = 0

#tlaková škála
def H(T, lnP, z):
	return np.exp(lnP)/(rho(T, lnP, z)*g(z))
#	return 1e3*(R*T)/(mu(T,lnP)*g(z))



#	integrace rovnic
#____________________________________________________________

def derivative(y, z):
	T, lnP = y
	r = eos.rhoofp(X, ztab, T*1e-6, np.exp(lnP)*1e-11, 1)	#rhoofp vraci hustotu v g/cm^3
	reteos = eos.esac(X, ztab, T*1e-6, r, 6, irad)
	P, E, S, cv, chiRho, chiT = reteos[:6]
#	return np.array([ T/(chiT*H(T, lnP, z)), 1/H(T,lnP,z) ])
	return np.array([ T*TDeriv(T, lnP, z)/H(T, lnP, z), 1/H(T,lnP,z) ])
#	return np.array([ (ext.T(z+100)-ext.T(z))/100., (np.log(ext.p(z+100))-np.log(ext.p(z)))/100. ])
#	return np.array([ (ext.T(z+100)-ext.T(z))/100., 1/H(T, lnP, z) ])


#počáteční podmínky
P0 = (2./3.)*g(0)/opac(ext.T(0), np.log(ext.p(0)),0)
P0 = ext.p(0)
#Teff = 4000
Teff = ext.T(0)
yinit = np.array([ Teff, np.log(P0) ])

#integrace pomocí isoda pack
z = np.linspace(0., 12.5e6, 500)
y = integrate.odeint(derivative, yinit, z, mxstep=500000)
T = y[:,0]
Tfunc = interpolate.interp1d(z,T)
lnP = y[:,1]
pe = np.exp(lnP)

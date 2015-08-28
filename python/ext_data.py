# -*- coding: utf-8 -*-
import numpy as np
from interp_wrapper import interp
from scipy import interpolate
from scipy.interpolate import UnivariateSpline
execfile('/home/tom/Documents/skola/bakalarka/python/const.py')

z = np.genfromtxt("/home/tom/Documents/skola/bakalarka/model_S2.dat", usecols = (0))
arr_size = z.size
r = np.empty(arr_size)
r[0] = 0
for i in range (1,arr_size):
    r[i] = r[i-1]+(z[arr_size-i]-z[arr_size-i-1])

temperature = np.genfromtxt("/home/tom/Documents/skola/bakalarka/model_S2.dat", usecols = (1))
pressure = np.genfromtxt("/home/tom/Documents/skola/bakalarka/model_S2.dat", usecols = (2))
rho = np.genfromtxt("/home/tom/Documents/skola/bakalarka/model_S2.dat", usecols = (3))
opac = np.genfromtxt("/home/tom/Documents/skola/bakalarka/model_S2.dat", usecols = (4))
hscale = np.genfromtxt("/home/tom/Documents/skola/bakalarka/model_S2.dat", usecols = (10))
adGrad = np.genfromtxt("/home/tom/Documents/skola/bakalarka/model_S2.dat", usecols = (11))
delta = np.genfromtxt("/home/tom/Documents/skola/bakalarka/model_S2.dat", usecols = (12))

mass = np.empty(arr_size)
mass[0] = 0

for i in range (1,arr_size):
    mass[i] = mass[i-1]+4*np.pi*rho[arr_size-i-1]*(r[i]*r[i])*(r[i]-r[i-1])

gravity = np.zeros(arr_size)

for i in range (0,arr_size-1):
    gravity[i] = G*mass[arr_size-i-1]/(r[arr_size-i-1]*r[arr_size-i-1])

Z = np.linspace(0,12.5e6,500)

T = UnivariateSpline(z,temperature,s=0)
p = UnivariateSpline(z,pressure,s=0)
g = UnivariateSpline(z,gravity,s=0)
r = UnivariateSpline(z,rho,s=0)
o = UnivariateSpline(z,opac,s=0)
H = UnivariateSpline(z,hscale,s=0)
adG = UnivariateSpline(z,adGrad,s=0)
adG_T = UnivariateSpline(temperature,adGrad,s=0)
d = UnivariateSpline(z,delta,s=0)
M = UnivariateSpline(z,mass[::-1],s=0)

#tabulka kombinací http://webs.wichita.edu/physics/opacity/ a http://cdsweb.u-strasbg.fr/topbase/TheOP.html

op = np.genfromtxt('/home/tom/Documents/skola/bakalarka/python/opacity.dat')

logR = np.linspace(-8., 1., 19)
logT = np.linspace(3., 6., 61)

f = interpolate.RectBivariateSpline(logT, logR, op, kx=2, ky=2)

def opac(T, rho, z):
	T6 = T*1e-6
	R = (rho*1e-3)/T6**3
#	if R>10:
#		print "R  =  "
#		print R
	log_opac = f(np.log10(T),np.log10(R))
	return 10**log_opac[0,0]*0.1
#	return o(z)

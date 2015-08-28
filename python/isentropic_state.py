import os.path
import sys
sys.path.append('/home/tom/Documents/skola/bakalarka/python/')
sys.path.append('/home/tom/Documents/skola/bakalarka/python/EOS/')

import numpy as np
from scipy.optimize import bisect
from OPAL.EOS import gas_state
import EOS
import EOS_id
from ext_data import g
import ext_data as ext
from const import *
import pdb

def derivative(z, rho, cp, cv, S):
	gamma = cp/cv
	ret = (g(z)/gamma)*rho**(2-gamma)*np.exp(-gamma*S/cp)
	pdb.set_trace()
	return ret

z = np.linspace(0., 12.5e6, 2000)
dz = np.diff(z)
T = np.empty(z.size)
P = np.empty(z.size)
rho = np.empty(z.size)
mu = np.empty(z.size)
cp = np.empty(z.size)
cv = np.empty(z.size)

T[0] = 3500.
S0 = EOS.S(ext.T(12.5e6), ext.p(12.5e6), 12.5e6)
P[0] = bisect(lambda P: EOS.S(T[0],P,0)-S0,1e1,ext.p(0))
S0 = S0*1e2
state = gas_state(T[0], P[0])
cp[0] = state[1]
cv[0] = state[2]
mu[0] = EOS_id.mu(T[0], P[0])
rho[0] = EOS_id.rho(T[0], P[0], z[0])

model = [z[0], T[0], P[0]]
for i in range(1,z.size):
	T[i] = P[i-1]*mu[i-1]*1e-3/(rho[i-1]*R)
	state = gas_state(T[i], P[i-1])
	cp[i] = state[1]
	cv[i] = state[2]
	rho[i] = rho[i]+dz[i-1]*derivative(z[i], rho[i-1], cp[i], cv[i], S0)
	mu[i] = EOS_id.mu(T[i], P[i-1])
	P[i] = T[i]*rho[i]*R/(1e-3*mu[i])
	model.append([z[i], T[i], P[i]])

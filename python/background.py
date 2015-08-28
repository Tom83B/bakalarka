# -*- coding: utf-8 -*-
import numpy as np
from scipy import integrate
from scipy.integrate import ode
from scipy import interpolate
from scipy.optimize import bisect

import .ext_data as ext
from .ext_data import g
from .const import *
from .derivatives import derivative, logPderivative, TDeriv 
from .OPAL import EOS

irad = 0


#	integrace rovnic
#____________________________________________________________


#=======počáteční podmínky=======
S0 = S(ext.T(12.5e6), ext.p(12.5e6), 12.5e6)
top = 0.e6
bottom = 12.5e6
start = top
end = bottom if start==top else top
T0 = ext.T(start)
P0 = ext.p(start)
#P0 = bisect(lambda P: S(T0,P,0)-S0,1e1,ext.p(0))
y0 = np.array([ start, np.log(T0) ])

#=========integrace================
model = []
f = lambda logP, y: logPderivative(y,logP)
r = ode(f).set_integrator('dopri5')
r.set_initial_value(y0, np.log(P0))
dt = 0.01
#while r.successful() and r.t<bottom:
while r.successful() and r.y[0]<bottom:
	z = r.y[0]
	print z
	T = np.exp(r.y[1])
	P = np.exp(r.t)
	new_state = np.concatenate(([z], [T], [P], EOS.gas_state(T,P)))
	new_state = np.append(new_state, TDeriv(T,P,r.t))
	model.append(new_state.tolist())
	r.integrate(r.t+dt)

model = np.asarray(model)

#   uložení do souboru
#_______________________________________________________________

#y[:,1]=np.exp(y[:,1])
#data = np.concatenate((np.array([z]).T, y), axis=1)
#data = np.concatenate((np.array([z]).T, np.array([ext.T(z)]).T, np.array([ext.p(z)]).T), axis=1)
np.savetxt('background.dat', model)

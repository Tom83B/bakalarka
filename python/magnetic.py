# -*- coding: utf-8 -*-
import os.path
import sys
sys.path.append('/home/tom/Documents/skola/bakalarka/python/')

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from const import phi

import loadbcg as bcg
from const import Btop, phi

def mag_iteration(y, P, eps):
	d_main = np.empty(500)
	d_sub = np.zeros(500)
	d_sup = np.zeros(500)
	b = np.zeros(500)

	h = 25000	#grid spacing

	Bz0 = np.sqrt(8e-7*np.pi*(bcg.P[499]-P[499]))
	b[0] = np.sqrt(Btop)
	b[499] = np.sqrt(Bz0)

	for i in range(1,499): d_main[i] = y[i]**3+eps*((phi/np.pi)*y[i]/h**2+y[i]**3)
	for i in range(0,498): d_sub[i] = -(eps*phi*y[i+1])/(2*np.pi*h*h)
	for i in range(2,500): d_sup[i] = -(eps*phi*y[i-1])/(2*np.pi*h*h)
	for i in range(1,499): b[i] = y[i]**4+(eps*8e-7*np.pi*(bcg.P[i]-P[i]))
	d_main[0] = 1
	d_main[499] = 1

	data = [d_sub, d_main, d_sup]
	diags = [-1,0,1]
	A = sparse.spdiags(data, diags, 500, 500, format='csc')

	return spsolve(A,b)

def error(v, w):
	return np.linalg.norm(v-w)/np.linalg.norm(v)

def comp_mag(y, P):
	err = 1
	y_new = y
	while err>1e-2:
		err_diff = 2
		i = 0
		err_new = err
		
		while err_new >= err:
			y_new = mag_iteration(y, P, 2**(-i))
			err_new = error(y, y_new)
			i = i+1
			print i

		y = y_new
		err = err_new
	return y_new
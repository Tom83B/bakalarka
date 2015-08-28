# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

z = np.genfromtxt("/home/tom/Documents/škola/bakalářka/model_S.dat", usecols = (0))
arr_size = z.size
r = np.empty(arr_size)
r[0] = 0
for i in range (1,arr_size):
    r[i] = r[i-1]+(z[arr_size-i]-z[arr_size-i-1])
    print "kurva"
    if i==arr_size-1:
        print "ahoj"

T_ext = np.genfromtxt("/home/tom/Documents/škola/bakalářka/model_S.dat", usecols = (1))
p_ext = np.genfromtxt("/home/tom/Documents/škola/bakalářka/model_S.dat", usecols = (2))
rho_ext = np.genfromtxt("/home/tom/Documents/škola/bakalářka/model_S.dat", usecols = (3))

mass = np.empty(arr_size)
mass[0] = 0

for i in range (1,arr_size+1):
    mass[i] = mass[i-1]+4*np.pi*rho_ext[arr_size-i-1]*(r[i]*r[i])*(r[i]-r[i-1])

g = np.zeros(arr_size)

for i in range (0,arr_size-1):
    g[i] = G*mass[arr_size-i-1]/(r[arr_size-i-1]*r[arr_size-i-1])


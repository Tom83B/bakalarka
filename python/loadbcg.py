# -*- coding: utf-8 -*-
import numpy as np

z = np.genfromtxt("/home/tom/Documents/skola/bakalarka/python/background.dat", usecols = (0))
T = np.genfromtxt("/home/tom/Documents/skola/bakalarka/python/background.dat", usecols = (1))
P = np.genfromtxt("/home/tom/Documents/skola/bakalarka/python/background.dat", usecols = (2))
rho = np.genfromtxt("/home/tom/Documents/skola/bakalarka/python/background.dat", usecols = (3))
cp = np.genfromtxt("/home/tom/Documents/skola/bakalarka/python/background.dat", usecols = (4))
cv = np.genfromtxt("/home/tom/Documents/skola/bakalarka/python/background.dat", usecols = (5))
nabla_ad = np.genfromtxt("/home/tom/Documents/skola/bakalarka/python/background.dat", usecols = (6))
delta = np.genfromtxt("/home/tom/Documents/skola/bakalarka/python/background.dat", usecols = (7))
S = np.genfromtxt("/home/tom/Documents/skola/bakalarka/python/background.dat", usecols = (8))
cp2 = np.genfromtxt("/home/tom/Documents/skola/bakalarka/python/background.dat", usecols = (9))

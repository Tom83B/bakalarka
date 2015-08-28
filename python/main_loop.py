from derivatives import Frad, Fconv
from EOS import rho, cp
from const import dt

def disconnection_check(y):
	B = y*y
	return False

while disconnection_check(y) == False:
	for i in range (0, z.size):
		FconvArr[i] = Fconv(T[i], P[i], z[i])
		FradArr[i] = Frad(T[i], P[i], z[i])
		rhoArr[i] = rho(T[i], P[i], z[i])
		cpArr[i] = cp(T[i], P[i])

	F = FconvArr+FradArr
	dFdz = np.diff(F)/h
	dFdz = np.append(dFdz, dFdz[dFdz.size-1])
	
	dTdt = -dFdz/(cpArr*rhoArr)
	T = T+dt*dTdt
	P0 = adjust_baseP(y, z, T, rho(T[499], P[499], z[499]))
	P = integrate.odeint(lambda P, z: PDeriv(T_func(z), P, z), P0, z_inv)

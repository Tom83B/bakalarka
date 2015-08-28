def get_flows():
	FradArr = np.empty(z.size)
	FconvArr = np.empty(z.size)
	extFradArr = np.empty(z.size)
	extFconvArr = np.empty(z.size)
	
	for i in range (0, z.size):
		extFconvArr[i] = Fconv(ext.T(z[i]), np.log(ext.p(z[i])), z[i])
		extFradArr[i] = Frad(ext.T(z[i]), np.log(ext.p(z[i])), z[i])
		FconvArr[i] = Fconv(y[i,0], y[i,1], z[i])
		FradArr[i] = Frad(y[i,0], y[i,1], z[i])
	
	return np.array([ FradArr, FconvArr, extFradArr, extFconvArr ])

def get_opacity():
	oArr = np.empty(z.size)
	extoArr = np.empty(z.size)

	for i in range (0, z.size):
		oArr[i] = opac(T[i], lnP[i], z[i])
		extoArr[i] = opac(ext.T(z[i]), np.log(ext.p(z[i])), z[i])
	
	return np.array([ oArr, extoArr ])

def get_rho():
	rhoArr = np.empty(z.size)
	for i in range (0, z.size):
		rhoArr[i] = rho(T[i], lnP[i], z[i])
	
	return rhoArr

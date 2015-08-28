
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

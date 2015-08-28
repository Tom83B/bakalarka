def get_model(T0, P0, top=0., bottom=12.5e6, isentropic=False, alpha=1.9904568, dt=0.01)
	model = []
	def f(logP, y):
		state = EOS.gas_state(T,P)
		z = y[0]
		T = np.exp(y[1])
		rho = state[0]
		adG = state[3]
		H = np.exp(logP)/(rho*g(z))
		return [ H, TDeriv(T, P, state) if isentropic==False else adG ]
	f = lambda logP, y: [H]
	r = ode(f).set_integrator('dopri5')
	r.set_initial_value(y0, np.log(P0))
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
	return model

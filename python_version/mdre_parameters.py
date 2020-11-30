import numpy as np


class Parameters:
	
	# constructor
	def __init__(self):
		pass
	
	sc        = 1.0
	K_inj     = 0.0							# optical injection strength
	E0        = 10.0
	omega_inj = 1.0
	
	
	# laser parameters
	
	a_L              = 15					# number of active QW-layers
	f_act            = 0.5					# fraction of resonant (active) QDs
	f_inact          = 1.0 - f_act			# fraction of off-resonant (inactive) QDs
	g_GS             = 230.0				# ground state gain coefficient in ns⁻¹
	N_QD             = 1.0					# 2D InAs QD density per QW layer
	R_w_loss         = 0.54					# QW loss rate in ns⁻¹
	T_e_eq           = 300.0				# quasi-equilibrium electron temperature in K
	T_h_eq           = 300.0				# quasi-equilibrium hole temperature in K
	W_GS             = 0.44					# coefficient for spontaneous recombination losses of ground state in ns⁻¹
	W_ES             = 0.55					# coefficient for spontaneous recombination losses of excited state in ns⁻¹
	
	kappa            = 50.0					# cavity loss rate in ns⁻¹
	delta_omega_ES   = 125.0				# change of instantaneous frequency (excited state) in ns⁻¹
	delta_omega_e_QW = 11.3					# change of instantaneous frequency (electrons) in ns⁻¹
	delta_omega_h_QW = 5.5					# change of instantaneous frequency (holes) in ns⁻¹
	eps_ES_e         = -14.0				# electron excited state confinement energy in meV
	eps_GS_e         = -64.0				# electron ground state confinement energy in meV
	eps_GS_h         = -35.0				# hole ground state confinement energy in meV
	eps_ES_h         = -15.0				# hole excited state confinement energy in meV
	
	
	# constants
	
	e      = 1.6021766208e-19				# elementary charge in C
	hbar   = 6.582119514e-13				# Planck constant in meVs
	kB     = 8.6173303e-2					# Boltzmann constant in meVK⁻¹
	m0     = 9.10938356e-31					# electron mass in kg
	me     = 0.043 * m0						# effective electron mass in kg
	mh     = 0.45 * m0						# effective hole mass in kg
	
	# 2D density of states
	
	Dens_e = me / (np.pi * np.power(hbar, 2))
	Dens_h = mh / (np.pi * np.power(hbar, 2))
	
	# fit parameters for charge-carrier scattering processes
	
	A_GS_e = 18.5							# in ns⁻¹
	B_GS_e = 1.9
	A_GS_h = 10.5							# in ns⁻¹
	B_GS_h = 5.3
	
	A_ES_e = 48.3							# in ns⁻¹
	B_ES_e = 0.48
	A_ES_h = 21.4							# in ns⁻¹
	B_ES_h = 1.8
	
	C_e    = 1014.0							# in ns⁻¹
	D_e    = 1.4
	C_h    = 2272.0							# in ns⁻¹
	D_h    = 2.3
	

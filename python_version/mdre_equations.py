import numpy as np
#~ import sympy as sy

import mdre_parameters
import mdre_processes as pr

p = mdre_parameters.Parameters()


# right hand sides of MDRE model:

# A real electric field amplitude
def eq0 (A, phi, rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e, J=17.0):
	return ((pr.g(rho_GS_e_act, rho_GS_h_act) - p.kappa) * A)


# phi electric field phase
def eq1 (A, phi, rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e, J=17.0):
	return (-pr.delta_omega(rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e))


# rho_GS_e_act electron occupation probability of active QD in ground state
def eq2 (A, phi, rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e, J=17.0):
	return (-pr.g(rho_GS_e_act, rho_GS_h_act) * np.power(A, 2) / (p.f_act * p.a_L * p.N_QD) - pr.R_sp_GS_act(rho_GS_e_act, rho_GS_h_act) + pr.S_GS_e_cap_act(rho_GS_e_act, w_e) + pr.S_e_rel_act(rho_GS_e_act, rho_ES_e, w_e))


# rho_GS_h_act hole occupation probability of active QD in ground state
def eq3 (A, phi, rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e, J=17.0):
	return (-pr.g(rho_GS_e_act, rho_GS_h_act) * np.power(A, 2) / (p.f_act * p.a_L * p.N_QD) - pr.R_sp_GS_act(rho_GS_e_act, rho_GS_h_act) + pr.S_GS_h_cap_act(rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e) + pr.S_h_rel_act(rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e))


# rho_GS_e_inact electron occupation probability of inactive QD in ground state
def eq4 (A, phi, rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e, J=17.0):
	return (-pr.R_sp_GS_inact(rho_GS_e_inact, rho_GS_h_inact) + pr.S_GS_e_cap_inact(rho_GS_e_inact, w_e) + pr.S_e_rel_inact(rho_GS_e_inact, rho_ES_e, w_e))


# rho_GS_h_inact hole occupation probability of inactive QD in ground state
def eq5 (A, phi, rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e, J=17.0):
	return (-pr.R_sp_GS_inact(rho_GS_e_inact, rho_GS_h_inact) + pr.S_GS_h_cap_inact(rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e) + pr.S_h_rel_inact(rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e))


# rho_ES_e electron occupation probability of excited state
def eq6 (A, phi, rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e, J=17.0):
	return (-pr.R_sp_ES(rho_ES_e, rho_ES_h) + pr.S_ES_e_cap(rho_ES_e, w_e) - 0.5 * (p.f_act * pr.S_e_rel_act(rho_GS_e_act, rho_ES_e, w_e) + p.f_inact * pr.S_e_rel_inact(rho_GS_e_inact, rho_ES_e, w_e)))


# rho_ES_h hole occupation probability of excited state
def eq7 (A, phi, rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e, J=17.0):
	return (-pr.R_sp_ES(rho_ES_e, rho_ES_h) + pr.S_ES_h_cap(rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e) - 0.5 * (p.f_act * pr.S_h_rel_act(rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e) + p.f_inact * pr.S_h_rel_inact(rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e)))


# w_e 2D electron density
def eq8 (A, phi, rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e, J=17.0):
	return (-2.0 * p.N_QD * (p.f_act * pr.S_GS_e_cap_act(rho_GS_e_act, w_e) + p.f_inact * pr.S_GS_e_cap_inact(rho_GS_e_inact, w_e)) - 4.0 * p.N_QD * pr.S_ES_e_cap(rho_ES_e, w_e) - pr.r_loss_w(rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e) + J)


# w_h 2D hole density
def eq9 (A, phi, rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e, J=17.0):
	return (-2.0 * p.N_QD * (p.f_act * pr.S_GS_h_cap_act(rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e) + p.f_inact * pr.S_GS_h_cap_inact(rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e)) - 4.0 * p.N_QD * pr.S_ES_h_cap(rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e) - pr.r_loss_w(rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e) + J)



# returns right hand side of MDRE model
def mdre_RHS (x0, J = 60.0):
	
	A              = x0[0]
	phi            = x0[1]
	rho_GS_e_act   = x0[2]
	rho_GS_h_act   = x0[3]
	rho_GS_e_inact = x0[4]
	rho_GS_h_inact = x0[5]
	rho_ES_e       = x0[6]
	rho_ES_h       = x0[7]
	w_e            = x0[8]
	
	dA              = eq0(A, phi, rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e)
	dphi            = eq1(A, phi, rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e)
	drho_GS_e_act   = eq2(A, phi, rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e)
	drho_GS_h_act   = eq3(A, phi, rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e)
	drho_GS_e_inact = eq4(A, phi, rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e)
	drho_GS_h_inact = eq5(A, phi, rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e)
	drho_ES_e       = eq6(A, phi, rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e)
	drho_ES_h       = eq7(A, phi, rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e)
	dw_e            = eq8(A, phi, rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e, J)
	
	return (np.array([dA, dphi, drho_GS_e_act, drho_GS_h_act, drho_GS_e_inact, drho_GS_h_inact, drho_ES_e, drho_ES_h, dw_e]))



# returns symbolic state vector of MDRE model with reduced system dimension
#~ def mdre_state_vector ():
    
    #~ A, phi                         = sy.symbols("A, phi")
    #~ rho_GS_e_act, rho_GS_h_act     = sy.symbols("rho_GS_e_act, rho_GS_h_act")
    #~ rho_GS_e_inact, rho_GS_h_inact = sy.symbols("rho_GS_e_inact, rho_GS_h_inact")
    #~ rho_ES_e, rho_ES_h             = sy.symbols("rho_ES_e, rho_ES_h")
    #~ w_e                            = sy.symbols("w_e")
    
    #~ return (sy.Matrix([A, phi, rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e]))



# returns symbolic right hand side of MDRE model with reduced system dimension
# w_h is replaced with w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact))
#~ def mdre_sym_RHS ():
    
    #~ # symbolic dynamic variables
    #~ A, phi                         = sy.symbols("A, phi")
    #~ rho_GS_e_act, rho_GS_h_act     = sy.symbols("rho_GS_e_act, rho_GS_h_act")
    #~ rho_GS_e_inact, rho_GS_h_inact = sy.symbols("rho_GS_e_inact, rho_GS_h_inact")
    #~ rho_ES_e, rho_ES_h             = sy.symbols("rho_ES_e, rho_ES_h")
    #~ w_e                            = sy.symbols("w_e")
    
    #~ # symbolic right hand side of MDRE model with reduced system dimension
    #~ dmdre = sy.Matrix([(p.g_GS * (rho_GS_e_act + rho_GS_h_act - 1.0) - p.kappa) * A, 
                       #~ -p.delta_omega_ES * (rho_ES_e + rho_ES_h) + p.delta_omega_e_QW * w_e + p.delta_omega_h_QW * (w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact))), 
                       #~ -p.g_GS * (rho_GS_e_act + rho_GS_h_act - 1.0) * A * A / (p.f_act * p.a_L * p.N_QD) - p.W_GS * rho_GS_e_act * rho_GS_h_act + p.A_GS_e * w_e * w_e / (p.B_GS_e + w_e) * (1.0 - rho_GS_e_act) - p.A_GS_e * w_e * w_e / (p.B_GS_e + w_e) * sy.exp((p.eps_GS_e - p.kB * p.T_e_eq * sy.log(sy.exp((w_e * p.e * 1e12) / (p.Dens_e * p.kB * p.T_e_eq)) - 1.0)) / (p.kB * p.T_e_eq)) * rho_GS_e_act + p.C_e * w_e / (p.D_e + w_e) * rho_ES_e * (1.0 - rho_GS_e_act) - p.C_e * w_e / (p.D_e + w_e) * sy.exp((p.eps_GS_e - p.eps_ES_e) / (p.kB * p.T_e_eq)) * (1.0 - rho_ES_e) * rho_GS_e_act, 
                       #~ -p.g_GS * (rho_GS_e_act + rho_GS_h_act - 1.0) * A * A / (p.f_act * p.a_L * p.N_QD) - p.W_GS * rho_GS_e_act * rho_GS_h_act + (p.A_GS_h * (w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact))) * (w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact))) / (p.B_GS_h + w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact)))) * (1.0 - rho_GS_h_act) - ((p.A_GS_h * (w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact))) * (w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact))) / (p.B_GS_h + w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact)))) * sy.exp((p.eps_GS_h - (p.kB * p.T_h_eq * sy.log(sy.exp(((w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact))) * p.e * 1e12) / (p.Dens_h * p.kB * p.T_h_eq)) - 1.0))) / (p.kB * p.T_h_eq))) * rho_GS_h_act + (p.C_h * (w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact))) / (p.D_h + (w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact))))) * rho_ES_h * (1.0 - rho_GS_h_act) - ((p.C_h * (w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact))) / (p.D_h + (w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact))))) * sy.exp((p.eps_GS_h - p.eps_ES_h) / (p.kB * p.T_h_eq))) * (1.0 - rho_ES_h) * rho_GS_h_act, 
                       #~ -p.W_GS * rho_GS_e_inact * rho_GS_h_inact + p.A_GS_e * w_e * w_e / (p.B_GS_e + w_e) * (1.0 - rho_GS_e_inact) - p.A_GS_e * w_e * w_e / (p.B_GS_e + w_e) * sy.exp((p.eps_GS_e - p.kB * p.T_e_eq * sy.log(sy.exp((w_e * p.e * 1e12) / (p.Dens_e * p.kB * p.T_e_eq)) - 1.0)) / (p.kB * p.T_e_eq)) * rho_GS_e_inact + p.C_e * w_e / (p.D_e + w_e) * rho_ES_e * (1.0 - rho_GS_e_inact) - p.C_e * w_e / (p.D_e + w_e) * sy.exp((p.eps_GS_e - p.eps_ES_e) / (p.kB * p.T_e_eq)) * (1.0 - rho_ES_e) * rho_GS_e_inact, 
                       #~ -p.W_GS * rho_GS_e_inact * rho_GS_h_inact + (p.A_GS_h * (w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact))) * (w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact))) / (p.B_GS_h + w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact)))) * (1.0 - rho_GS_h_inact) - ((p.A_GS_h * (w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact))) * (w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact))) / (p.B_GS_h + w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact)))) * sy.exp((p.eps_GS_h - p.kB * p.T_h_eq * sy.log(sy.exp(((w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact))) * p.e * 1e12) / (p.Dens_h * p.kB * p.T_h_eq)) - 1.0)) / (p.kB * p.T_h_eq))) * rho_GS_h_inact + (p.C_h * (w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact))) / (p.D_h + (w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact))))) * rho_ES_h * (1.0 - rho_GS_h_inact) - ((p.C_h * (w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact))) / (p.D_h + (w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact))))) * sy.exp((p.eps_GS_h - p.eps_ES_h) / (p.kB * p.T_h_eq))) * (1.0 - rho_ES_h) * rho_GS_h_inact, 
                       #~ -p.W_ES * rho_ES_e * rho_ES_h + p.A_ES_e * w_e * w_e / (p.B_ES_e + w_e) * (1.0 - rho_ES_e) - p.A_ES_e * w_e * w_e / (p.B_ES_e + w_e) * sy.exp((p.eps_ES_e - p.kB * p.T_e_eq * sy.log(sy.exp((w_e * p.e * 1e12) / (p.Dens_e * p.kB * p.T_e_eq)) - 1.0)) / (p.kB * p.T_e_eq)) * rho_ES_e - 0.5 * (p.f_act * (p.C_e * w_e / (p.D_e + w_e) * rho_ES_e * (1.0 - rho_GS_e_act) - p.C_e * w_e / (p.D_e + w_e) * sy.exp((p.eps_GS_e - p.eps_ES_e) / (p.kB * p.T_e_eq)) * (1.0 - rho_ES_e) * rho_GS_e_act) + p.f_inact * (p.C_e * w_e / (p.D_e + w_e) * rho_ES_e * (1.0 - rho_GS_e_inact) - p.C_e * w_e / (p.D_e + w_e) * sy.exp((p.eps_GS_e - p.eps_ES_e) / (p.kB * p.T_e_eq)) * (1.0 - rho_ES_e) * rho_GS_e_inact)), 
                       #~ -p.W_ES * rho_ES_e * rho_ES_h + (p.A_ES_h * (w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact))) * (w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact))) / (p.B_ES_h + (w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact))))) * (1.0 - rho_ES_h) - ((p.A_ES_h * (w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact))) * (w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact))) / (p.B_ES_h + (w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact))))) * sy.exp((p.eps_ES_h - p.kB * p.T_h_eq * sy.log(sy.exp(((w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact))) * p.e * 1e12) / (p.Dens_h * p.kB * p.T_h_eq)) - 1.0)) / (p.kB * p.T_h_eq))) * rho_ES_h - 0.5 * (p.f_act * ((p.C_h * (w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact))) / (p.D_h + (w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact))))) * rho_ES_h * (1.0 - rho_GS_h_act) - (p.C_h * (w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact))) / (p.D_h + (w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact))))) * sy.exp((p.eps_GS_h - p.eps_ES_h) / (p.kB * p.T_h_eq)) * (1.0 - rho_ES_h) * rho_GS_h_act) + p.f_inact * ((p.C_h * (w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact))) / (p.D_h + (w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact))))) * rho_ES_h * (1.0 - rho_GS_h_inact) - ((p.C_h * (w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact))) / (p.D_h + (w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact))))) * sy.exp((p.eps_GS_h - p.eps_ES_h) / (p.kB * p.T_h_eq))) * (1.0 - rho_ES_h) * rho_GS_h_inact)), 
                       #~ -2.0 * p.N_QD * (p.f_act * (p.A_GS_e * w_e * w_e / (p.B_GS_e + w_e) * (1.0 - rho_GS_e_act) - p.A_GS_e * w_e * w_e / (p.B_GS_e + w_e) * sy.exp((p.eps_GS_e - p.kB * p.T_e_eq * sy.log(sy.exp((w_e * p.e * 1e12) / (p.Dens_e * p.kB * p.T_e_eq)) - 1.0)) / (p.kB * p.T_e_eq)) * rho_GS_e_act) + p.f_inact * (p.A_GS_e * w_e * w_e / (p.B_GS_e + w_e) * (1.0 - rho_GS_e_inact) - p.A_GS_e * w_e * w_e / (p.B_GS_e + w_e) * sy.exp((p.eps_GS_e - p.kB * p.T_e_eq * sy.log(sy.exp((w_e * p.e * 1e12) / (p.Dens_e * p.kB * p.T_e_eq)) - 1.0)) / (p.kB * p.T_e_eq)) * rho_GS_e_inact)) - 4.0 * p.N_QD * (p.A_ES_e * w_e * w_e / (p.B_ES_e + w_e) * (1.0 - rho_ES_e) - p.A_ES_e * w_e * w_e / (p.B_ES_e + w_e) * sy.exp((p.eps_ES_e - p.kB * p.T_e_eq * sy.log(sy.exp((w_e * p.e * 1e12) / (p.Dens_e * p.kB * p.T_e_eq)) - 1.0)) / (p.kB * p.T_e_eq)) * rho_ES_e) - p.R_w_loss * w_e * (w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact)))])
    
    #~ return (dmdre)

























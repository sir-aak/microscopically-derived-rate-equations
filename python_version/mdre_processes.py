import numpy as np
import mdre_parameters

p = mdre_parameters.Parameters()

# due to charge conservation and to achieve a non-singular Jacobian matrix
# the charge carrier density w_h is replaced as follows:
# w_h = (w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact)))


# amplitude gain in ns⁻¹
def g (rho_GS_e_act, rho_GS_h_act):
	return (p.g_GS * (rho_GS_e_act + rho_GS_h_act - 1.0))


# change of instantaneous frequency in ns⁻¹
def delta_omega (rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e):
	return (p.delta_omega_ES * (rho_ES_e + rho_ES_h) + p.delta_omega_e_QW * w_e + p.delta_omega_h_QW * (w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact))))


# spontaneous emission term
def dE_sp (E, rho_GS_e_act, rho_GS_h_act):
	return (p.beta * p.a_L * p.N_QD * p.f_act * p.W_GS * rho_GS_e_act * rho_GS_h_act * E / np.power(np.absolute(E), 2))


# *********************************************************************


# electron quasi Fermi level
# p.e * 1e12 is a compensation factor that ensures the usage of the correct physical units
def E_F_e_eq (w_e):
	return (p.kB * p.T_e_eq * np.log(np.exp((w_e * p.e * 1e12) / (p.Dens_e * p.kB * p.T_e_eq)) - 1.0))


# hole quasi Fermi level
# p.e * 1e12 is a compensation factor that ensures the usage of the correct physical units
def E_F_h_eq (rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e):
	return (p.kB * p.T_h_eq * np.log(np.exp(((w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact))) * p.e * 1e12) / (p.Dens_h * p.kB * p.T_h_eq)) - 1.0))


# *********************************************************************


# spontaneous recombination losses

def R_sp_GS_act (rho_GS_e_act, rho_GS_h_act):
	return (p.W_GS * rho_GS_e_act * rho_GS_h_act)


def R_sp_GS_inact (rho_GS_e_inact, rho_GS_h_inact):
	return (p.W_GS * rho_GS_e_inact * rho_GS_h_inact)


def R_sp_ES (rho_ES_e, rho_ES_h):
	return (p.W_ES * rho_ES_e * rho_ES_h)


# effective QW loss rate
def r_loss_w (rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e):
	return (p.R_w_loss * w_e * (w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact))))


# *********************************************************************


# nonlinear scattering rates:

# in-scattering rates derived from full microscopically calculated rates

# capture in-rates
def S_GS_e_cap_in (w_e):
	return (p.A_GS_e * np.power(w_e, 2) / (p.B_GS_e + w_e))


def S_GS_h_cap_in (rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e):
	return (p.A_GS_h * np.power(w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact)), 2) / (p.B_GS_h + w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact))))


def S_ES_e_cap_in (w_e):
	return (p.A_ES_e * np.power(w_e, 2) / (p.B_ES_e + w_e))


def S_ES_h_cap_in (rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e):
	return (p.A_ES_h * np.power(w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact)), 2) / (p.B_ES_h + (w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact)))))


# relaxation in-rates
def S_e_rel_in (w_e):
	return (p.C_e * w_e / (p.D_e + w_e))


def S_h_rel_in (rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e):
	return (p.C_h * (w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact))) / (p.D_h + (w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact)))))


# *********************************************************************


# out-scattering rates given by detailed balance relationships

# capture out-rates
def S_GS_e_cap_out (w_e):
	
	if w_e == 0.0:
		return (0.0)
	
	else:
		return (S_GS_e_cap_in(w_e) * np.exp((p.eps_GS_e - E_F_e_eq(w_e)) / (p.kB * p.T_e_eq)))


def S_GS_h_cap_out (rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e):
	
	if ((w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact))) == 0.0):
		return (0.0)
	
	else:
		return (S_GS_h_cap_in(rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e) * np.exp((p.eps_GS_h - E_F_h_eq(rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e)) / (p.kB * p.T_h_eq)))


def S_ES_e_cap_out (w_e):
	
	if (w_e == 0.0):
		return (0.0)
	
	else:
		return (S_ES_e_cap_in(w_e) * np.exp((p.eps_ES_e - E_F_e_eq(w_e)) / (p.kB * p.T_e_eq)))


def S_ES_h_cap_out (rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e):
	
	if ((w_e + 4.0 * p.N_QD * (rho_ES_e - rho_ES_h) + 2.0 * p.N_QD * (p.f_act * (rho_GS_e_act - rho_GS_h_act) + p.f_inact * (rho_GS_e_inact - rho_GS_h_inact))) == 0.0):
		return (0.0)
	
	else:
		return (S_ES_h_cap_in(rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e) * np.exp((p.eps_ES_h - E_F_h_eq(rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e)) / (p.kB * p.T_h_eq)))


# relaxation out-rates
def S_e_rel_out (w_e):
	return (S_e_rel_in(w_e) * np.exp((p.eps_GS_e - p.eps_ES_e) / (p.kB * p.T_e_eq)))


def S_h_rel_out (rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e):
	return (S_h_rel_in(rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e) * np.exp((p.eps_GS_h - p.eps_ES_h) / (p.kB * p.T_h_eq)))



# *********************************************************************


# scattering processes between different charge-carrier states:

# capture rates
def S_GS_e_cap_act (rho_GS_e_act, w_e):
	return (S_GS_e_cap_in(w_e) * (1.0 - rho_GS_e_act) - S_GS_e_cap_out(w_e) * rho_GS_e_act)


def S_GS_h_cap_act (rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e):
	return (S_GS_h_cap_in(rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e) * (1.0 - rho_GS_h_act) - S_GS_h_cap_out(rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e) * rho_GS_h_act)


def S_GS_e_cap_inact (rho_GS_e_inact, w_e):
	return (S_GS_e_cap_in(w_e) * (1.0 - rho_GS_e_inact) - S_GS_e_cap_out(w_e) * rho_GS_e_inact)


def S_GS_h_cap_inact (rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e):
	return (S_GS_h_cap_in(rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e) * (1.0 - rho_GS_h_inact) - S_GS_h_cap_out(rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e) * rho_GS_h_inact)


def S_ES_e_cap (rho_ES_e, w_e):
	return (S_ES_e_cap_in(w_e) * (1.0 - rho_ES_e) - S_ES_e_cap_out(w_e) * rho_ES_e)


def S_ES_h_cap (rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e):
	return (S_ES_h_cap_in(rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e) * (1.0 - rho_ES_h) - S_ES_h_cap_out(rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e) * rho_ES_h)


# relaxation rates
def S_e_rel_act (rho_GS_e_act, rho_ES_e, w_e):
	return (S_e_rel_in(w_e) * rho_ES_e * (1.0 - rho_GS_e_act) - S_e_rel_out(w_e) * (1.0 - rho_ES_e) * rho_GS_e_act)


def S_h_rel_act (rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e):
	return (S_h_rel_in(rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e) * rho_ES_h * (1.0 - rho_GS_h_act) - S_h_rel_out(rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e) * (1.0 - rho_ES_h) * rho_GS_h_act)


def S_e_rel_inact (rho_GS_e_inact, rho_ES_e, w_e):
	return (S_e_rel_in(w_e) * rho_ES_e * (1.0 - rho_GS_e_inact) - S_e_rel_out(w_e) * (1.0 - rho_ES_e) * rho_GS_e_inact)


def S_h_rel_inact (rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e):
	return (S_h_rel_in(rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e) * rho_ES_h * (1.0 - rho_GS_h_inact) - S_h_rel_out(rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e) * (1.0 - rho_ES_h) * rho_GS_h_inact)

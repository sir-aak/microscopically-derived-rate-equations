import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


# plots linear fit function to DnuInj_crit(beta) data for injection
def plotLinearFitInjection ():
	
	# parameters
	Ttrans = 100
	Teval  = 250
	K_inj  = 0.1
	
	# load data for linear fit
	data = np.loadtxt("DnuInj_crit(beta)_Ttrans=100_Teval=250_linear.txt")
	
	# sorts data for first column
	data   = data[data[:, 0].argsort()]
	beta   = data[:, 0]
	DnuInj = data[:, 1]
	
	def linearFit (x, m, n):
		return (m * x + n)
	
	# initial guess parameters for linear fit:
	m = 0.033
	n = -1.2591
	initial_guess = [m, n]
	
	# fitting function parameters to data
	popt, pcov = curve_fit(linearFit, beta, DnuInj, p0=initial_guess)
	
	m = np.around(popt[0], 5)
	n = np.around(popt[1], 5)
	
	m_err = np.around(np.sqrt(pcov[0, 0]), 6)
	n_err = np.around(np.sqrt(pcov[1, 1]), 6)
	
	plt.figure(figsize=(4.6, 3.6))
	# plt.suptitle("MDRE: SNIPER bifurcation parameter", fontsize=14)
	# plt.title(r"$K_{inj}=$" + str(K_inj) + r", $T_{trans}=$" + str(Ttrans) + r"$\,$ns, $T_{eval}=$" + str(Teval)+ r"$\,$ns")
	plt.subplots_adjust(top=0.98, bottom=0.15, left=0.23, right=0.98)
	plt.semilogx(beta, DnuInj, color="black", marker=".", markersize=5, linewidth=0)
	plt.semilogx(beta, linearFit(beta, *popt), color="red")
	# ~ plt.semilogx(beta, linearFit(beta, *initial_guess), color="gray")
	plt.xlabel(r"spontaneous emission coefficient $\beta$")
	plt.ylabel(r"critical master-slave detuning $\Delta\nu_{inj,crit}$ / GHz")
	plt.text(2e-7, -1.229, r"fit function: $\Delta\nu_{inj,crit}(\beta) = m\beta + n$")
	plt.text(2e-7, -1.2315, r"fitted parameters:")
	plt.text(2e-7, -1.234, r"$m=($" + str(m) + r" $\pm$ " + str(m_err) + r") GHz")
	plt.text(2.5e-7, -1.2365, r"$n=($" + str(n) + r" $\pm$ " + str(n_err) + r") GHz")
	plt.grid(color="lightgray")
	plt.show()


# plots power fit function to DnuInj_crit(beta) data for injection
def plotPowerFitInjection ():
	
	# parameters
	Ttrans = 100
	Teval  = 250
	K_inj  = 0.1
	
	# load data for power fit
	data = np.loadtxt("DnuInj_crit(beta)_Ttrans=100_Teval=250_power.txt")
	
	# sorts data for first column
	data = data[data[:, 0].argsort()]
	beta = data[:, 0]
	
	# offset for injection at K_inj = 0.1, corresponding DnuInjcrit at beta=0
	DnuInj = data[:, 1] + 1.25867405274989 #1.2586740
	
	def powerFit (x, a, b):
		return (a * x ** b)
	
	# initial guess parameters for power fit:
	a = 0.033
	b = 1.0
	initial_guess = [a, b]
	
	# fitting function parameters to data
	popt, pcov = curve_fit(powerFit, beta, DnuInj, p0=initial_guess)
	
	a = np.around(popt[0], 5)
	b = np.around(popt[1], 4)
	
	a_err = np.around(np.sqrt(pcov[0, 0]), 5)
	b_err = np.around(np.sqrt(pcov[1, 1]), 4)
	
	
	plt.figure(figsize=(4.6, 3.6))
	# plt.suptitle("MDRE: SNIPER bifurcation parameter", fontsize=14)
	# plt.title(r"$K_{inj}=$" + str(K_inj) + r", $T_{trans}=$" + str(Ttrans) + r"$\,$ns, $T_{eval}=$" + str(Teval)+ r"$\,$ns")
	plt.subplots_adjust(top=0.98, bottom=0.15, left=0.23, right=0.98)
	plt.loglog(beta, DnuInj, color="black", marker=".", markersize=5, linewidth=0)
	plt.loglog(beta, powerFit(beta, *popt), color="red")
	# ~ plt.loglog(beta, powerFit(beta, *initial_guess), color="gray")
	plt.xlabel(r"spontaneous emission coefficient $\beta$")
	plt.ylabel(r"critical master-slave detuning difference" + "\n" + r"$\Delta\nu_{inj,crit}(\beta) - \Delta\nu_{inj,crit}(\beta=0)$ / GHz")
	plt.text(1.2e-5, 1.5e-2, r"fit function: $\Delta\nu_{inj,crit}(\beta) = a\,\beta^b$")
	plt.text(1.2e-5, 6e-3, r"fitted parameters:")
	plt.text(1.2e-5, 2.6e-3, r"$a=($" + str(a) + r" $\pm$ " + str(a_err) + r") GHz")
	plt.text(1.2e-5, 1.1e-3, r"$b=$" + str(b) + r" $\pm$ " + str(b_err))
	plt.grid(color="lightgray")
	plt.show()


# plots power fit function to DnuInj_crit(beta) data 
# for injection with delayed feedback
def plotPowerFitInjectionFeedback ():
	
	# parameters
	Ttrans = 1000
	Teval  = 250
	K_fb   = 0.05
	tau    = 0.08
	phi    = np.pi / 2.0
	K_inj  = 0.1
	
	# load data for power fit
	data = np.loadtxt("DnuInj_crit(beta)_injection_power.txt")
	
	# sorts data for first column
	# ~ data = data[data[:, 0].argsort()]
	beta = data[:, 0]
	
	# offset for injection and feedback at K_inj = 0.1
	DnuInj = data[:, 1] + 1.25867405274989#1.8850417
	
	def powerFit (x, a, b):
		return (a * x ** b)
	
	# initial guess parameters for power fit:
	a = 0.033
	b = 1.0
	initial_guess = [a, b]
	
	# fitting function parameters to data
	popt, pcov = curve_fit(powerFit, beta, DnuInj, p0=initial_guess)
	
	a = np.around(popt[0], 5)
	a_err = np.around(np.sqrt(pcov[0, 0]), 5)
	
	b = np.around(popt[1], 4)
	b_err = np.around(np.sqrt(pcov[1, 1]), 4)
	
	textypos = 0.9
	
	plt.figure(figsize=(5.9, 4.6))
	ax = plt.gca()
	plt.rcParams.update({"font.size": 14.0})
	# plt.suptitle("MDRE: SNIPER bifurcation parameter", fontsize=14)
	# plt.title(r"$K_{inj}=$" + str(K_inj) + r", $T_{trans}=$" + str(Ttrans) + r"$\,$ns, $T_{eval}=$" + str(Teval)+ r"$\,$ns")
	plt.subplots_adjust(top=0.96, bottom=0.16, left=0.24, right=0.98)
	plt.loglog(beta, DnuInj, color="black", marker=".", markersize=5, linewidth=0)
	plt.loglog(beta, powerFit(beta, *popt), color="red")
	# ~ plt.loglog(beta, powerFit(beta, *initial_guess), color="gray")
	plt.xlabel(r"spontaneous emission coefficient $\beta$", fontsize=18.0)
	plt.ylabel(r"difference of critical detunings" + "\n" + r"$\Delta\nu_{inj,crit}(\beta) - \Delta\nu_{inj,crit}(\beta=0)$ / GHz", fontsize=18.0)
	plt.xticks([0.0000001, 0.00001, 0.001, 0.1], fontsize=18.0)
	plt.yticks([0.01, 0.0001, 0.000001, 0.00000001], fontsize=18.0)
	plt.text(0.06, textypos, r"fit function: $\Delta\nu_{inj,crit}(\beta) = a\,\beta^b$", transform = ax.transAxes)
	plt.text(0.06, textypos-0.075, r"fitted parameters:", transform = ax.transAxes)
	plt.text(0.06, textypos-2*0.075, r"$a=$" + str(a) + r" $\pm$ " + str(a_err), transform = ax.transAxes)
	plt.text(0.06, textypos-3*0.075, r"$b=$" + str(b) + r" $\pm$ " + str(b_err), transform = ax.transAxes)
	plt.grid(color="lightgray")
	plt.show()


# ~ plotLinearFitInjection()
# ~ plotPowerFitInjection()
plotPowerFitInjectionFeedback()


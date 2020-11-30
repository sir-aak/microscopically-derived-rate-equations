import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


# plots SNIPER bifurcation scan of MDRE model
def plotSNIPER ():
	
	data      = np.loadtxt("sniper.txt")
	mu        = data[:, 0]
	period    = data[:, 1]
	amplitude = data[:, 2]
	
	# parameter at which SNIPER bifurcation occurs
	DnuInjCrit = -1.2588
	
	def func (x, a, b):
		return (a * x**b)
	
	a = 0.8
	b = -0.5
	initial_guess = [a, b]
	
	popt, pcov = curve_fit(func, mu, period, p0 = initial_guess)
	
	a     = np.around(popt[0], 4)
	b     = np.around(popt[1], 4)
	a_err = np.around(np.sqrt(pcov[0, 0]), 4)
	b_err = np.around(np.sqrt(pcov[1, 1]), 4)
	
	# ~ fitted_parameters = "fitted parameters:\n" + "$a$ = " + str(np.around(popt[0], 3)) + "\n$b$ = " + str(np.around(popt[1], 3))
	
	plt.figure(figsize=(5.9, 5.9))
	plt.subplots_adjust(top=0.97, bottom=0.17, left=0.19, right=0.98)
	plt.rcParams.update({"font.size": 18})
	plt.loglog(mu, period, color="black", marker="o", markersize=5, label=r"$T(\mu)$")
	plt.loglog(mu, func(mu, *popt), color="red", label="fitted function")
	plt.xlabel(r"bifurcation parameter" + "\n" + r"$\mu = \Delta\nu_{inj,crit} - \Delta\nu_{inj}$ / GHz", size=18)
	plt.ylabel(r"period of limit cycle $T$ / ns", size=18)
	plt.ylim(1.0, 100.0)
	plt.text(0.005, 65.0, r"$\Delta\nu_{inj,crit} = $" + str(DnuInjCrit) + r"$\,$GHz")
	plt.text(0.005, 45.0, r"fit function: $f(\mu) = a\,\mu^b$")
	plt.text(0.005, 33.0, "fit parameters:")
	plt.text(0.005, 24.0, r"$a=$" + str(a) + r"$\pm$" + str(a_err))
	plt.text(0.005, 18.0, r"$b=$" + str(b) + r"$\pm$" + str(b_err))
	# ~ plt.grid(color="lightgray")
	plt.legend(loc="lower left")
	
	plt.figure(figsize=(5.9, 5.9))
	plt.subplots_adjust(top=0.97, bottom=0.17, left=0.19, right=0.98)
	plt.loglog(mu, amplitude, color="black", marker="o", markersize=5)
	plt.xlabel(r"bifurcation parameter" + "\n" + r"$\mu = \Delta\nu_{inj,crit} - \Delta\nu_{inj}$ / GHz", size=18)
	plt.ylabel(r"amplitude of limit cycle $|E|^2$", size=18)
	plt.ylim(1e-1, 1e1)
	# ~ plt.grid(color="lightgray")
	
	plt.show()


# plots Hopf bifurcation scan of MDRE model
def plotHopf ():
	
	data      = np.loadtxt("hopf.txt")
	mu        = data[:, 0]
	period    = data[:, 1]
	amplitude = data[:, 2]
	
	# parameter at which Hopf bifurcation occurs
	DnuInjCrit = 2.647
	
	plt.figure(figsize=(5.9, 5.9))
	plt.subplots_adjust(top=0.97, bottom=0.17, left=0.25, right=0.97)
	plt.rcParams.update({"font.size": 18})
	plt.loglog(mu, period, color="black", marker="o", markersize=5)
	plt.text(0.004, 0.5, r"$\Delta\nu_{inj,crit} = $" + str(DnuInjCrit) + r"$\,$GHz")
	plt.xlabel(r"bifurcation parameter" + "\n" + r"$\mu = \Delta\nu_{inj,crit} - \Delta\nu_{inj}$ / GHz")
	plt.ylabel(r"period of limit cycle $T$ / ns")
	plt.ylim(1e-1, 1.0)
	plt.grid(color="lightgray")
	
	plt.figure(figsize=(5.9, 5.9))
	plt.subplots_adjust(top=0.97, bottom=0.17, left=0.25, right=0.97)
	plt.plot(mu, amplitude, color="black", marker="o", markersize=5)
	plt.xlabel(r"bifurcation parameter" + "\n" + r"$\mu = \Delta\nu_{inj,crit} - \Delta\nu_{inj}$ / GHz", labelpad=10)
	plt.ylabel(r"amplitude of limit cycle $|E|^2$")
	plt.grid(color="lightgray")
	
	plt.show()


# ~ plotSNIPER()
plotHopf()


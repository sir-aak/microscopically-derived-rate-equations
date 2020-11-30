import numpy as np
import matplotlib.pyplot as plt


# plots free running laser frequency and electric field amplitude 
# as a function of beta for MDRE model
def plotFrequency ():
	
	data   = np.loadtxt("free_running_laser.txt")
	beta   = data[:, 0]
	omega0 = data[:, 1]
	E0     = data[:, 2]
	
	# ~ plt.figure(figsize=(4.6, 3.6))
	# ~ plt.subplots_adjust(top=0.98, bottom=0.15, left=0.17, right=0.98)
	# ~ plt.subplots_adjust(top=0.98, left=0.14, right=0.98)
	# ~ plt.suptitle("MDRE: frequency of free running laser", fontsize=14)
	# ~ plt.semilogx(beta, omega0, color="black")
	# ~ plt.xlabel(r"spontaneous emission coefficient $\beta$")
	# ~ plt.ylabel(r"frequency of free running laser $\omega_0$ / ns$^{-1}$")
	# ~ plt.grid(color="lightgray")
	# ~ plt.show()
	
	fig, ax1 = plt.subplots(1, 1)
	ax2      = ax1.twinx()
	fig.set_size_inches(5.9, 3.9)
	plt.rcParams.update({"font.size": 12.5})
	fig.subplots_adjust(top=0.98, bottom=0.16, left=0.19, right=0.82)
	ax1.semilogx(beta, omega0, color="black")
	ax2.semilogx(beta, E0, color="red")
	ax1.set_xlabel(r"spontaneous emission coefficient $\beta$", fontsize=12.5)
	ax1.set_ylabel(r"angular frequency" + "\n" 
	              + "of free running laser $\omega_0$ / ns$^{-1}$", fontsize=12.5)
	ax2.set_ylabel(r"electric field amplitude" + "\n" 
	              + "of free running laser $E_0$", color="red", fontsize=12.5)
	ax1.tick_params(axis="x", labelsize=12.5)
	ax1.tick_params(axis="y", labelsize=12.5)
	ax2.yaxis.label.set_color('red')
	ax2.tick_params(axis="y", colors="red", labelsize=12.5)
	# ~ ax.tick_params(axis="x", labelsize=12)
	# ~ ax2.set_xticklabels(fontsize=12.5)
	# ~ ax1.set_yticklabels([162.6, 162.8, 163.0, 163.2, 163.4, 163.6, 163.8], fontsize=12.5)
	plt.show()

plotFrequency()


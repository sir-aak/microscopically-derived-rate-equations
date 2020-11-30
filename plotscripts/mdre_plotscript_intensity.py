import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d


# plots intensity time series for MDRE model
def plotIntensity ():
	
	# index boundaries for time 3D plot
	nStart = 0
	nEnd   = 250000
	
	with open("time_series.txt", "r") as file:
		lines = file.readlines()
	
	time         = []
	intensity    = []
	# ~ rho_GS_e_act = []
	# ~ rho_GS_h_act = []
	E_real       = []
	E_imag       = []
	
	for line in lines:
		time.append(float((line.split('	')[0])))
		intensity.append(float((line.split('	')[1])))
		E_real.append(float((line.split('	')[2])))
		E_imag.append(float((line.split('	')[3])))
		# ~ rho_GS_e_act.append(float((line.split('	')[6])))
		# ~ rho_GS_h_act.append(float((line.split('	')[7])))
	
	time         = np.array(time)
	intensity    = np.array(intensity)
	E_real       = np.array(E_real)
	E_imag       = np.array(E_imag)
	# ~ rho_GS_e_act = np.array(rho_GS_e_act)
	# ~ rho_GS_h_act = np.array(rho_GS_h_act)
	
	# ~ inversion = rho_GS_e_act + rho_GS_h_act - 1.0
	
	plt.figure(figsize=(9, 6))
	plt.suptitle("MDRE: intensity time series", fontsize=14)
	plt.title(r"$K_{inj} = 0.05$, $\Delta\nu_{inj} = -5$, $K_{fb} = 0.085$, $\varphi = 2.5$, $\tau = 0.8\,$ns")
	plt.plot(time, intensity)
	# ~ plt.plot(time, inversion)
	#~ plt.xlim(1040.0, 1050.0)
	#~ plt.ylim(0.3, 2.1)
	plt.xlabel(r"time $t$ / ns")
	plt.ylabel(r"intensity $|E|^2$")
	plt.grid(color="lightgray")
	
	# ~ plt.figure()
	# ~ plt.rcParams.update({"font.size": 12})
	# ~ ax = plt.axes(projection='3d')
	# ~ ax.plot3D(time[nStart:nEnd], E_real[nStart:nEnd], E_imag[nStart:nEnd])
	# ~ ax.set_xlabel(r"time $t$ / ns")
	# ~ ax.set_ylabel(r"Re($E$)")
	# ~ ax.set_zlabel(r"Im($E$)")
	
	plt.show()


plotIntensity()


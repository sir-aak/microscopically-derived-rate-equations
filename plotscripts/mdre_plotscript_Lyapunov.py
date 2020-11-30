import numpy as np
import matplotlib.pyplot as plt


# plots intensity time series for MDRE model
def plotLyapunovExponent ():
	
	# parameters
	K_inj  = 0.2
	DnuInj = -4
	K_fb   = 0.05
	phi    = 0
	tau    = 0.8
	
	with open("time_series.txt", "r") as file:
		lines = file.readlines()
	
	time      = []
	intensity = []
	Lyapunov  = []
	
	for line in lines:
		time.append(float((line.split('	')[0])))
		intensity.append(float((line.split('	')[1])))
		Lyapunov.append(float((line.split('	')[22])))
	
	time      = np.array(time)
	intensity = np.array(intensity)
	Lyapunov  = np.array(Lyapunov)
	
	
	plt.figure(figsize=(9, 6))
	plt.suptitle("MDRE: Lyapunov exponent time series", fontsize=14)
	plt.title(r"$K_{inj}=$" + str(K_inj) + r", $\Delta\nu_{inj}=$" + str(DnuInj) + r", $K_{fb}=$" + str(K_fb) + r", $\varphi=$" + str(phi) + r", $\tau=$" + str(tau) + r"$\,$ns")
	plt.plot(time, Lyapunov)
	# ~ plt.xlim(5.0, 20.0)
	# ~ plt.ylim(np.min(Lyapunov[5000:20000]), np.max(Lyapunov[5000:20000]))
	# ~ plt.text(12.0, -90.0, "Lyapunov exponent = " + str(0.008))
	plt.xlabel(r"time $t$ / ns")
	plt.ylabel(r"intensity $|E|^2$")
	plt.grid(color="lightgray")
	
	plt.show()


plotLyapunovExponent()


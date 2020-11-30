import numpy as np
import matplotlib.pyplot as plt


# plots occupation time series for MDRE model
def plotOccupation ():
	
	with open("time_series.txt", "r") as file:
		lines = file.readlines()
	
	time          = []
	intensity     = []
	occupation_GS = []
	occupation_ES = []
	
	for line in lines:
		time.append(float((line.split('	')[0])))
		intensity.append(float((line.split('	')[1])))
		occupation_GS.append(float((line.split('	')[16])))
		occupation_ES.append(float((line.split('	')[17])))
	
	y_shift = 0.43
	
	time          = np.array(time)
	intensity     = np.array(intensity)
	occupation_GS = np.array(occupation_GS)
	occupation_ES = np.array(occupation_ES) + y_shift
	
	
	fig, ax1 = plt.subplots()
	fig.set_size_inches(10.0, 7.0)
	fig.suptitle("MDRE occupation time series", fontsize=14)
	
	ax1.set_title(r"@ $K=0.06$, $\Delta\nu_{inj}=0.95$ GHz")
	p1 = ax1.plot(time, intensity, color="red", alpha=0.5, label=r"intensity $|E|^2$")[0]
	ax1.set_ylabel(r"intensity $|E|^2$")
	ax1.set_xlabel(r"time $t$ / ns")
	
	ax2 = ax1.twinx()
	p2 = ax2.plot(time, occupation_GS, color="red", label=r"ground state $\rho_{GS,e}^{act} + \rho_{GS,h}^{act}$")[0]
	p3 = ax2.plot(time, occupation_ES, color="blue", label=r"excited state $\rho_{ES,e} + \rho_{ES,h}$ + " + str(y_shift))[0]
	ax2.set_ylabel(r"occupation $\rho_e + \rho_h$")
	
	plots = [p1, p2, p3]
	ax2.legend(plots, [p.get_label() for p in plots], loc="lower right", fancybox=True, framealpha=0.85)
	
	plt.show()


plotOccupation()

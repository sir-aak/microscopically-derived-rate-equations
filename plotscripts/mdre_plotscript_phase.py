import numpy as np
import matplotlib.pyplot as plt


# plots phase time series for MDRE model
def plotPhase ():
	
	with open("time_series.txt", "r") as file:
		lines = file.readlines()
	
	time  = []
	phase = []
	
	for line in lines:
		
		time.append(float((line.split('	')[0])))
		phase.append(float((line.split('	')[3])))
	
	time  = np.array(time)
	phase = np.array(phase)
	
	plt.figure()
	plt.suptitle("MDRE phase time series", fontsize=14)
	plt.plot(time, phase)
	plt.xlabel(r"time $t$ / ns")
	plt.ylabel(r"$\varphi(t)$")
	plt.grid(color="lightgray")
	plt.show()


plotPhase()

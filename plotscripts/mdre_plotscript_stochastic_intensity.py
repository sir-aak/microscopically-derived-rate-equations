import numpy as np
import matplotlib.pyplot as plt


# plots deterministic and stochastic intensity time series for MDRE model
def plotIntensity ():
	
	# get data from deterministic time series
	
	time      = []
	intensity = []
	
	# ~ with open("time_series.txt", "r") as file:
		# ~ lines = file.readlines()
	
	# ~ for line in lines:
		# ~ time.append(float((line.split('	')[0])))
		# ~ intensity.append(float((line.split('	')[1])))
	
	# ~ time      = np.array(time)
	# ~ intensity = np.array(intensity)
	
	
	# get data from stochastic time series
	
	time_stoch      = []
	intensity_stoch = []
	
	with open("time_series_stochastic.txt", "r") as file:
		lines = file.readlines()
	
	for line in lines:
		time_stoch.append(float((line.split('	')[0])))
		intensity_stoch.append(float((line.split('	')[1])))
	
	time_stoch      = np.array(time_stoch)
	intensity_stoch = np.array(intensity_stoch)
	
	
	plt.figure(figsize=(8, 7))
	plt.suptitle("MDRE intensity time series", fontsize=14)
	plt.title(r"$K=0.1$, $\Delta\nu_{inj}=0$, $\beta=10^{-2}$, d$t=10^{-3}\,$ns, method: Euler-Maruyama")
	plt.plot(time_stoch, intensity_stoch, color="blue", label="stochastic")
	# ~ plt.plot(time, intensity, color="red", label="deterministic")
	plt.hlines([0.65, 0.7], np.min(time_stoch), np.max(time_stoch))
	plt.xlabel(r"time $t$ / ns")
	plt.ylabel(r"intensity $|E|^2$")
	plt.legend(loc="lower right")
	plt.grid(color="lightgray")
	plt.show()


plotIntensity()


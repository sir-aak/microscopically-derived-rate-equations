import numpy as np
import matplotlib.pyplot as plt


# plots deterministic and stochastic intensity time series for MDRE model
def plotIntensity ():
	
	nStart = 0
	nEnd   = 305000
	
	# get data from deterministic time series
	
	time      = []
	intensity = []
	
	with open("time_series_beta=5e-2.txt", "r") as file:
		lines = file.readlines()
	
	for line in lines:
		time.append(float((line.split('	')[0])))
		intensity.append(float((line.split('	')[1])))
	
	time      = np.array(time)
	intensity = np.array(intensity)
	
	
	# get data from stochastic time series
	
	time_stoch      = []
	intensity_stoch = []
	
	with open("time_series_stochastic_beta=5e-2.txt", "r") as file:
		lines = file.readlines()
	
	for line in lines:
		time_stoch.append(float((line.split('	')[0])))
		intensity_stoch.append(float((line.split('	')[1])))
	
	time_stoch      = np.array(time_stoch)
	intensity_stoch = np.array(intensity_stoch)
	
	plt.figure(figsize=(5.9, 5.5))
	ax = plt.gca()
	# ~ plt.rcParams.update({"font.size": 18})
	plt.subplots_adjust(top=0.99, bottom=0.16, left=0.14, right=0.96)
	# ~ plt.suptitle("MDRE intensity time series", fontsize=14)
	# ~ plt.title(r"$K=0.1$, $\Delta\nu_{inj}=0$, $\beta=10^{-2}$, d$t=10^{-3}\,$ns, method: Euler-Maruyama")
	# ~ plt.plot(time_stoch[nStart:nEnd], intensity_stoch[nStart:nEnd], color="darkblue", label="stochastic")
	plt.plot(time_stoch[nStart:nEnd*10], intensity_stoch[nStart:nEnd*10], color="darkblue", label="stochastic")
	plt.plot(time[nStart:nEnd], intensity[nStart:nEnd], color="red", linewidth=2.0, label="deterministic")
	ax.set_aspect(time[nEnd] / 2.25)
	plt.xticks([0.0, 100.0, 200.0, 300.0], fontsize=27.0)
	plt.yticks([0.0, 1.0, 2.0], fontsize=27.0)
	plt.xlabel(r"time $t$ / ns", fontsize=27.0)
	plt.ylabel(r"intensity $|E|^2$", fontsize=27.0)
	plt.xlim(-5.0, 305.0)
	plt.ylim(0.0, 2.25)
	# ~ plt.ylim(-0.05, 2.78)
	# ~ plt.legend(loc="lower right")
	# ~ plt.grid(color="lightgray")
	plt.show()


plotIntensity()


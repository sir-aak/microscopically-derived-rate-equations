import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors


# generates 1D-histogram for intensity at a specified time value t for MDRE system
# t-value (float) at which densityPlot() is evaluated must exist in time_series_stochastic.txt
# bins (int) are the number of classes for the histogram
def densityPlot1D (t, bins):
	
	# get data from stochastic time series
	with open("time_series_stochastic.txt", "r") as file:
		lines  = file.readlines()
	
	time      = []
	intensity = []
	
	for line in lines:
		time.append(float((line.split('	')[0])))
		intensity.append(float((line.split('	')[1])))
	
	time      = np.array(time)
	intensity = np.array(intensity)
	
	numberOfTimeSeries = (time == time[0]).sum()
	
	
	plt.figure(figsize=(8, 6))
	plt.suptitle(r"histogram for " + str(numberOfTimeSeries) + " time series at $t = $" + str(t), fontsize=14)
	plt.title(r"$K=0.1$, $\Delta\nu_{inj}=-2.0$, $\beta=5\cdot10^{-2}$, d$t=10^{-3}\,$ns, method: strong order 1.0")
	plt.hist(intensity[np.where(time == t)], bins)
	plt.xlabel(r"intensity $|E|^2$")
	plt.ylabel("frequency")
	plt.show()



# generates 2D-histogram for intensity for MDRE system
# xbins and ybins (int) are the number of classes for the histogram in x- and y-direction
def densityPlot2D (xbins, ybins):
	
	# parameters
	beta   = 0.01
	K_inj  = 0
	DnuInj = 0
	K_fb   = 0.2
	phi    = 1.5708
	tau    = 0.8
	
	# get data from deterministic time series
	with open("histograms/K_inj=0_DnuInj=0/time_series.txt", "r") as file:
		lines = file.readlines()
	
	time      = []
	intensity = []
	
	for line in lines:
		time.append(float((line.split('	')[0])))
		intensity.append(float((line.split('	')[1])))
	
	# get data from stochastic time series
	with open("histograms/K_inj=0_DnuInj=0/time_series_stochastic.txt", "r") as file:
		lines = file.readlines()
	
	time_stoch      = []
	intensity_stoch = []
	
	for line in lines:
		time_stoch.append(float((line.split('	')[0])))
		intensity_stoch.append(float((line.split('	')[1])))
	
	# conversion of lists to numpy arrays
	time            = np.array(time)
	intensity       = np.array(intensity)
	time_stoch      = np.array(time_stoch)
	intensity_stoch = np.array(intensity_stoch)
	
	numberOfTimeSeries = np.sum(time_stoch == time_stoch[0])
	numberOfSteps      = time.size
	
	# computation of mean and standard deviation of stochastic intensity
	intensity_mean = np.zeros(numberOfSteps)
	intensity_std  = np.zeros(numberOfSteps)
	
	for i in range(numberOfSteps):
		intensity_mean[i] = np.mean(intensity_stoch[i::numberOfSteps])
		intensity_std[i]  = np.std(intensity_stoch[i::numberOfSteps])
	
	# computation of root mean square deviation between deterministic and averaged stochastic intensity
	rms = np.around(np.sqrt(np.sum(np.square(intensity - intensity_mean))), 5)
	
	# time average after transient time Ttrans in case of fixed point dynamics
	Ttrans                    = 5.0
	skipSteps                 = int(Ttrans / (time[1] - time[0]))
	fixed_point_deterministic = np.around(intensity[-1], 5)
	fixed_point_stochastic    = np.around(np.mean(intensity_mean[skipSteps:]), 5)
	fixed_point_std           = np.around(np.std(intensity_mean[skipSteps:]), 5)
	
	
	#plot data
	plt.figure(figsize=(5.9, 4.6))
	plt.subplots_adjust(top=0.98, bottom=0.11, left=0.09, right=0.95)
	plt.rcParams.update({"font.size": 9})
	# ~ plt.suptitle(r"histogram for " + str(numberOfTimeSeries) + " time series", fontsize=14)
	# ~ plt.title(r"$\beta=$" + str(beta) + ", $K_{inj}=$" + str(K_inj) + r", $\Delta\nu_{inj}=$" + str(DnuInj) + r"$\,$GHz, $K_{fb}=$" + str(K_fb) + r", $\varphi=\frac{\pi}{2}$" + r", $\tau=$" + str(tau) + r"$\,$ns")
	plt.hist2d(time_stoch, intensity_stoch, bins=(xbins, ybins), norm=mpl.colors.LogNorm(), cmap="binary", zorder=2)
	plt.plot(time, intensity, color="red", alpha=0.667, label="deterministic", zorder=4)
	plt.plot(time, intensity_mean, color="blue", alpha=0.667, label="stochastic mean", zorder=3)
	plt.plot(time, intensity_mean + intensity_std, linestyle="--", color="blue", alpha=0.5, label="standard deviation", zorder=3)
	plt.plot(time, intensity_mean - intensity_std, linestyle="--", color="blue", alpha=0.5, zorder=3)
	
	# transform = axes[0, 0].transAxes ???
	plt.text(1.25, 2.10, r"rms error $I_{det} - I_{stoch}$ = " + str(rms), transform = axes[0, 0].transAxes)
	plt.text(1.25, 2.02, r"$I^*_{det}$   = " + str(fixed_point_deterministic), transform = axes[0, 0].transAxes)
	plt.text(1.25, 1.94, r"$I^*_{stoch}$ = " + str(fixed_point_stochastic) + "$\pm$" + str(fixed_point_std))
	
	plt.xlabel(r"time $t$ / ns")
	plt.ylabel(r"intensity $|E|^2$")
	plt.ylim((0.0, 2.25))
	plt.legend(loc="lower right")
	plt.grid(color="lightgray")
	cb = plt.colorbar()
	cb.set_label("frequency")
	plt.show()


#~ densityPlot1D(0.75, 100)
densityPlot2D(20000, 200)


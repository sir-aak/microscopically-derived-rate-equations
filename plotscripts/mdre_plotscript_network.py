import numpy as np
import matplotlib.pyplot as plt


# parameters
K_fb   = 0.328
tau    = 0.08
phi    = 2.86
K_inj  = 0
DnuInj = 0



def plotTimeSeriesSingleLaser ():
    
    data = np.loadtxt("time_series.txt")
    time      = data[:, 0]
    intensity = data[:, 1]
    
    plt.figure(figsize=(7, 5))
    plt.suptitle("MDRE: single laser", fontsize=14)
    plt.title(r"d$t=10^{-3}\,$ns, $\tau = 4\,$ns, $K = 0$, $\eta = 1$, $\varphi = 0$, method: Runge-Kutta 4")
    plt.plot(time, intensity)
    plt.xlabel(r"time $t$ / ns")
    plt.ylabel(r"intensity $|E|^2$")
    plt.grid(color="lightgray")
    plt.show()


def plotTimeSeriesSingleLaserComparison ():
    
    Ttrans = 5.0
    beta   = 0.01
    
    # get data from deterministic time series
    
    time_det      = []
    intensity_det = []
    
    with open("time_series.txt", "r") as file:
        Lines = file.readlines()
    
    for Line in Lines:
        time_det.append(float((Line.split('	')[0])))
        intensity_det.append(float((Line.split('	')[1])))
    
    time_det      = np.array(time_det)
    intensity_det = np.array(intensity_det)
    
    
    # get data from first time series
    
    time_1      = []
    intensity_1 = []
    
    with open("time_series_stochastic_dt=1e-3.txt", "r") as file:
        Lines = file.readlines()
    
    for Line in Lines:
        time_1.append(float((Line.split('	')[0])))
        intensity_1.append(float((Line.split('	')[1])))
    
    time_1      = np.array(time_1)
    intensity_1 = np.array(intensity_1)
    
    # computing mean and standard deviation of intensity after transient time Ttrans
    skipSteps_1      = int(Ttrans / (time_1[1] - time_1[0]))
    intensity_mean_1 = np.around(np.mean(intensity_1[skipSteps_1:]), 4)
    intensity_std_1  = np.around(np.std(intensity_1[skipSteps_1:]), 4)
    
    
    # get data from second time series
    
    time_2      = []
    intensity_2 = []
    
    with open("time_series_stochastic_dt=1e-4.txt", "r") as file:
        Lines = file.readlines()
    
    for Line in Lines:
        time_2.append(float((Line.split('	')[0])))
        intensity_2.append(float((Line.split('	')[1])))
    
    time_2      = np.array(time_2)
    intensity_2 = np.array(intensity_2)
    
    # computing mean and standard deviation of intensity after transient time Ttrans
    skipSteps_2      = int(Ttrans / (time_2[1] - time_2[0]))
    intensity_mean_2 = np.around(np.mean(intensity_2[skipSteps_2:]), 4)
    intensity_std_2  = np.around(np.std(intensity_2[skipSteps_2:]), 4)
    
    
    # get data from third time series
    
    time_3      = []
    intensity_3 = []
    
    with open("time_series_stochastic_dt=1e-5.txt", "r") as file:
        Lines = file.readlines()
    
    for Line in Lines:
        time_3.append(float((Line.split('	')[0])))
        intensity_3.append(float((Line.split('	')[1])))
    
    time_3      = np.array(time_3)
    intensity_3 = np.array(intensity_3)
    
    # computing mean and standard deviation of intensity after transient time Ttrans
    skipSteps_3      = int(Ttrans / (time_3[1] - time_3[0]))
    intensity_mean_3 = np.around(np.mean(intensity_3[skipSteps_3:]), 4)
    intensity_std_3  = np.around(np.std(intensity_3[skipSteps_3:]), 4)
    
    
    plt.figure(figsize=(9, 6))
    plt.suptitle("MDRE: complex gaussian noise variance check", fontsize=14)
    # plt.title(r"$K = 0$, $\Delta\nu_{inj} = 0\,$GHz, $\eta = 0.05$, $\varphi = 0$, $\tau = 0.08\,$ns")
    plt.title(r"var(Re($E$)) = $\frac{1}{2}$, var(Im($E$)) = $\frac{1}{2}$, $\beta=$" + str(beta) + ", $T_{trans}=$" + str(Ttrans))
    plt.plot(time_1, intensity_1, color="purple", alpha=0.75, label=r"d$t=10^{-3}\,$ns")
    plt.plot(time_2, intensity_2, color="blue", alpha=0.75, label=r"d$t=10^{-4}\,$ns")
    plt.plot(time_3, intensity_3, color="darkgreen", alpha=0.75, label=r"d$t=10^{-5}\,$ns")
    plt.plot(time_det, intensity_det, color="red", label="deterministic")
    plt.text(14.8, 1.28, r"$I^*_{det}$=" + str(np.around(intensity_det[-1], 4)))
    plt.text(13.0, 1.26, r"d$t=10^{-3}: \quad I^*=$" + str(intensity_mean_1) + r"$\pm$" + str(intensity_std_1))
    plt.text(13.0, 1.24, r"d$t=10^{-4}: \quad I^*=$" + str(intensity_mean_2) + r"$\pm$" + str(intensity_std_2))
    plt.text(13.0, 1.22, r"d$t=10^{-5}: \quad I^*=$" + str(intensity_mean_3) + r"$\pm$" + str(intensity_std_3))
    plt.xlabel(r"time $t$ / ns")
    plt.ylabel(r"intensity $|E|^2$")
    plt.xlim(5.0, 20.0)
    plt.ylim(0.9, 1.30)
    plt.legend(loc="lower right")
    plt.grid(color="lightgray")
    
    plt.show()


def plotTimeSeriesLaserNetwork ():
    
    data = np.loadtxt("time_series.txt")
    
    time = data[:, 0]
    I0   = data[:, 1]
    I1   = data[:, 2]
    I2   = data[:, 3]
    I3   = data[:, 4]
    
    plt.figure(figsize=(9, 6))
    # ~ plt.suptitle("MDRE: four node all-to-all network", fontsize=14)
    # ~ plt.title(r"$K_{fb}=$" + str(K_fb) + r", $\tau=$" + str(tau) + r"$\,$ns, $\varphi=$" + str(phi) + ", $K_{inj}=$" + str(K_inj) + r", $\Delta\nu_{inj}=$" + str(DnuInj) + r"$\,$GHz")
    plt.plot(time, I0, label="laser 0")
    plt.plot(time, I1, label="laser 1")
    plt.plot(time, I2, label="laser 2")
    plt.plot(time, I3, label="laser 3")
    plt.xlabel(r"time $t$ / ns")
    plt.ylabel(r"intensity of $i$-th laser in network $|E_i|^2$")
    # ~ plt.xlim(np.min(time) - 20.0, np.max(time) + 20.0)
    # ~ plt.ylim(0.0, np.max(I0) + 0.05)
    plt.legend(loc="upper right")
    plt.grid(color="lightgray")
    
    # ~ I0 = I0[150000:]
    # ~ I1 = I1[150000:]
    # ~ I2 = I2[150000:]
    # ~ I3 = I3[150000:]
    
    # ~ deviation01 = np.sqrt(np.mean(np.square(I0 - I1)))
    # ~ deviation02 = np.sqrt(np.mean(np.square(I0 - I2)))
    # ~ deviation03 = np.sqrt(np.mean(np.square(I0 - I3)))
    # ~ deviation12 = np.sqrt(np.mean(np.square(I1 - I2)))
    # ~ deviation13 = np.sqrt(np.mean(np.square(I1 - I3)))
    # ~ deviation23 = np.sqrt(np.mean(np.square(I2 - I3)))
    
    # ~ plt.text(197.0, 2.25, "rms-deviation $(150 - 200)\,$ns:")
    # ~ plt.text(197.0, 2.15, r"$d_{01}$ = " + str(deviation01))
    # ~ plt.text(197.0, 2.05, r"$d_{02}$ = " + str(deviation02))
    # ~ plt.text(197.0, 1.95, r"$d_{03}$ = " + str(deviation03))
    # ~ plt.text(197.0, 1.85, r"$d_{12}$ = " + str(deviation12))
    # ~ plt.text(197.0, 1.75, r"$d_{13}$ = " + str(deviation13))
    # ~ plt.text(197.0, 1.65, r"$d_{23}$ = " + str(deviation23))
    
    # ~ print("d01 = " + str(deviation01))
    # ~ print("d02 = " + str(deviation02))
    # ~ print("d03 = " + str(deviation03))
    # ~ print("d12 = " + str(deviation12))
    # ~ print("d13 = " + str(deviation13))
    # ~ print("d23 = " + str(deviation23))
    
    plt.show()


def plotTimeSeriesStochasticLaserNetwork ():
    
    data = np.loadtxt("time_series_stochastic.txt")
    
    time = data[:, 0]
    I0   = data[:, 1]
    I1   = data[:, 2]
    I2   = data[:, 3]
    I3   = data[:, 4]
    
    plt.figure(figsize=(9, 6))
    plt.suptitle("MDRE: four node all-to-all network", fontsize=14)
    # ~ plt.title(r"d$t=10^{-3}\,$ns, $\tau = 4\,$ns, $\beta = 10^{-2}$, $K = 0$, $\eta = 1$, $\varphi = 0$, method: Euler-Maruyama")
    plt.plot(time, I0, label="laser 0", alpha=0.85)
    plt.plot(time, I1, label="laser 1", alpha=0.85)
    plt.plot(time, I2, label="laser 2", alpha=0.85)
    plt.plot(time, I3, label="laser 3", alpha=0.85)
    plt.xlabel(r"time $t$ / ns")
    plt.ylabel(r"intensity of $i$-th laser in network $|E_i|^2$")
    plt.legend(loc="lower right")
    plt.grid(color="lightgray")
    
    plt.show()


# ~ plotTimeSeriesSingleLaser()
# ~ plotTimeSeriesSingleLaserComparison()
# ~ plotTimeSeriesLaserNetwork()
plotTimeSeriesStochasticLaserNetwork()


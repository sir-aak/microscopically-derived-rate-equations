import numpy as np
import matplotlib.pyplot as plt


def plotFeedbackBifurcationLineScan ():
	
	Ttrans = 1000
	Teval  = 50
	phi    = 2.35
	tau    = 0.08
	
	data    = np.loadtxt("feedback_bifurcation_line_scan_phi=3.141593_tau=0.800000_upsweep.txt")
	K_fb    = data[:, 0]
	extrema = data[:, 1]
	
	plt.figure(figsize=(5.9, 4.1))
	plt.subplots_adjust(top=0.98, bottom=0.12, left=0.11, right=0.98)
	# ~ plt.suptitle("MDRE: feedback bifurcation diagram", fontsize=14)
	# ~ plt.title(r"$\varphi=\pi$" + r", $\tau=$" + str(tau) + r"$\,$ns, $T_{trans}=$" + str(Ttrans) + r"$\,$ns, $T_{eval}=$" + str(Teval) + "$\,$ns")
	plt.scatter(K_fb, extrema, color="black", marker='o', s=0.30, zorder=2)
	plt.xlim(-0.005, 0.505)
	plt.ylim(-0.05, np.max(extrema) + 0.05)
	plt.xlabel(r"feedback coupling strength $K_{fb}$")
	plt.ylabel(r"intensity of unique extrema $|E_{extrema}|^2$")
	plt.grid(color="lightgray")
	plt.show()


def plotFeedbackBifurcationSweeps ():
	
	Ttrans = 1000
	Teval  = 50
	phi    = 2.35
	tau    = 0.8
	
	data       = np.loadtxt("feedback_bifurcation_line_scan_phi=3.141593_tau=0.800000_upsweep.txt")
	K_fb_up    = data[:, 0]
	extrema_up = data[:, 1]
	
	data         = np.loadtxt("feedback_bifurcation_line_scan_phi=3.141593_tau=0.800000_downsweep.txt")
	K_fb_down    = data[:, 0]
	extrema_down = data[:, 1]
	
	plt.figure(figsize=(5.9, 4.6))
	plt.subplots_adjust(top=0.97, bottom=0.15, left=0.12, right=0.97)
	plt.rcParams.update({"font.size": 18})
	# plt.suptitle("MDRE: feedback bifurcation diagram", fontsize=14)
	# plt.title(r"$\varphi=\pi$" + r", $\tau=$" + str(tau) + r"$\,$ns, $T_{trans}=$" + str(Ttrans) + r"$\,$ns, $T_{eval}=$" + str(Teval) + "$\,$ns")
	plt.scatter(K_fb_up, extrema_up, color="royalblue", marker='o', alpha=0.8, s=0.3, zorder=2, label="upsweep")
	plt.scatter(K_fb_down, extrema_down, color="red", marker='o', s=0.3, alpha=0.4, zorder=2, label="downsweep")
	plt.text(0.28, 3.6, r"$\tau=$" + str(tau) + r"$\,$ns")
	plt.text(0.28, 3.2, r"$\varphi=\pi$")
	plt.xlim(-0.005, 0.505)
	plt.ylim(-0.05, 4.05)
	plt.xticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5], ["0", "0.1", "0.2", "0.3", "0.4", "0.5"])
	plt.yticks([0.0, 1.0, 2.0, 3.0, 4.0])
	plt.xlabel(r"feedback coupling strength $K_{fb}$")
	plt.ylabel(r"intensity of unique extrema $|E|^2$")
	legend = plt.legend(loc="upper left")
	
	for lh in legend.legendHandles: 
		lh.set_alpha(1.0)
		lh._sizes = [10]
	
	plt.grid(color="lightgray")
	plt.show()


# ~ plotFeedbackBifurcationLineScan()
plotFeedbackBifurcationSweeps()


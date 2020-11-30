import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


# plots 2D injection diagrams for MDRE model: single laser
def plotInjectionSingleLaser ():
	
	# simulation parameters
	Ttrans = 1000
	Teval  = 50
	K_fb   = 0.05
	phi    = 0.0
	tau    = 0.8
	
	K_inj_line = 0.3
	
	# parameters for 257 x 385 map
	KinjStart   = 0.0
	KinjEnd     = 0.5
	KinjStep    = 0.001953125
	DnuInjStart = -6.0
	DnuInjEnd   = 6.0
	DnuInjStep  = 0.03125
	
	# initialize data and state matrix
	DnuInjSteps = int((DnuInjEnd - DnuInjStart) / DnuInjStep) + 1
	KinjSteps   = int((KinjEnd - KinjStart) / KinjStep) + 1
	
	cwIntensityData = np.zeros((DnuInjSteps, KinjSteps))
	maximaData      = np.zeros((DnuInjSteps, KinjSteps))
	periodData      = np.zeros((DnuInjSteps, KinjSteps))
	LyapunovData    = np.zeros((DnuInjSteps, KinjSteps))
	
	
	# read data and fill data matrix
	
	i = 0
	
	for K_inj in np.arange(KinjStart, KinjEnd + 1e-4, KinjStep):
		
		cwIntensity = []
		maxima      = []
		period      = []
		Lyapunov    = []
		
		with open("single_Laser/tau=8e-1/K_fb=0,05_phi=1,5pi_tau=8e-1/injection_line_scan_mdre_K_inj={:f}".format(K_inj) + ".txt", "r") as file:
			lines = file.readlines()
		
		for line in lines:
			cwIntensity.append(float((line.split('	')[2])))
			maxima.append(int((line.split('	')[3])))
			period.append(float((line.split('	')[5])))
			Lyapunov.append(float((line.split('	')[6])))
		
		cwIntensity = np.array(cwIntensity)
		maxima      = np.array(maxima)
		period      = np.array(period)
		Lyapunov    = np.array(Lyapunov)
		
		cwIntensityData[:, i] = cwIntensity
		maximaData[:, i]      = maxima
		periodData[:, i]      = period
		LyapunovData[:, i]    = Lyapunov
		
		i += 1
	
	
	# general title and labels
	suptitle    = "MDRE - single laser with delayed feedback: injection - "
	title       = r"$K_{fb}=$" + str(K_fb) + r", $\varphi=\frac{3\pi}{2}$" + r", $\tau=$" + str(tau) + r"$\,$ns, $T_{trans}=$" + str(Ttrans) + r"$\,$ns, $T_{eval}=$" + str(Teval) + r"$\,$ns"
	xlabel      = r"injection strength $K_{inj}$"
	ylabel      = r"master-slave detuning $\Delta\nu_{inj}$ / GHz"		#$\longrightarrow$
	xticks      = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
	xticksNames = ["0", "0.1", "0.2", "0.3", "0.4", "0.5"]
	
	# bound maxima
	maximaData[np.where(maximaData > 4)] = 5
	
	# define custom colors for ..
	
	# .. maxima plot
	mylightyellow = "#fff9bd"
	mydarkyellow  = "#fecb67"
	myorange      = "#fd8a3b"
	myred         = "#e41c1d"
	mykarmin      = "#be0126"
	
	# .. chaotic regions in period plot
	mygray        = "#bbbbbb"
	
	
	# define custom color map, bounds and ticks
	
	mycolormap     = mpl.colors.ListedColormap([mylightyellow, mydarkyellow, myorange, myred, mykarmin])
	colorBarTicks  = np.arange(5)
	colorBarBounds = np.arange(6) - 0.5
	
	mycolormap2    = mpl.colors.ListedColormap([mygray])
	
	extent = np.array([KinjStart, KinjEnd, DnuInjStart, DnuInjEnd])
	aspect = (KinjEnd - KinjStart) / (DnuInjEnd - DnuInjStart)
	
	
	# big magic is happening here
	mappable = mpl.cm.ScalarMappable(cmap=mycolormap)
	mappable.set_array(np.arange(5))
	
	
	# masking the arrays
	
	maximaDataBool  = maximaData == 0
	maximaData      = np.ma.masked_where(maximaDataBool == True, maximaData)
	cwIntensityData = np.ma.masked_where(maximaDataBool == False, cwIntensityData)
	
	maximaDataBool  = maximaData == 5
	chaoticRegion   = np.ones(periodData.shape)
	periodData      = np.ma.masked_where(maximaDataBool == True, periodData)
	chaoticRegion   = np.ma.masked_where(maximaDataBool == False, chaoticRegion)
	
	LyapunovBoolNeg = LyapunovData <= 0.0
	LyapunovDataPos = np.ma.masked_where(LyapunovBoolNeg == True, LyapunovData)
	LyapunovDataNeg = np.ma.masked_where(LyapunovBoolNeg == False, LyapunovData)
	
	'''
	# maxima diagram ###################################################
	
	fig, ax = plt.subplots()
	fig.set_size_inches(5.9, 5.0)
	fig.subplots_adjust(top=1.07, bottom=0.04, left=0.15, right=0.79)
	# ~ plt.suptitle(suptitle + "unique maxima")
	# ~ plt.title(title)
	maxima = plt.imshow(maximaData, origin="lower", extent=extent, cmap=mycolormap)
	cwInt  = plt.imshow(cwIntensityData, origin="lower", extent=extent, cmap="Blues_r")
	# ~ plt.plot([K_inj_line, K_inj_line], [-6.0, 6.0], color="black")
	plt.xlabel(xlabel, fontsize=18)
	plt.ylabel(ylabel, fontsize=18)
	plt.xticks(xticks, xticksNames, fontsize=18)
	plt.yticks(fontsize=18)
	ax.set_aspect(aspect)
	
	ax2 = fig.add_axes([0.82, 0.575, 0.03, 0.385])
	ax3 = fig.add_axes([0.82, 0.15, 0.03, 0.385])
	
	cb_int = fig.colorbar(cwInt, ax2)
	cb_int.ax.tick_params(labelsize=18)
	cb_int.set_label(r"cw intensity", fontsize=18, labelpad=20.0)  # -0.5 for s = 0.02, 6.5 for s = 1, 
	
	cb_max = fig.colorbar(mappable, ax3, boundaries=colorBarBounds, ticks=colorBarTicks)
	cb_max.ax.set_yticklabels(["1", "2", "3", "4", ">4"])
	cb_max.ax.tick_params(labelsize=18)
	cb_max.set_label("unique maxima", fontsize=18, labelpad=8.0) # 13.0 for s = 0.02, 8.5 for s = 1
	
	# ~ ax.scatter(0.1, -1.2588, color="black", s=20)
	# ~ ax.scatter(0.3, 2.647, color="red", s=20, zorder=3)
	# ~ ax.scatter(0.3, 2.8, marker="x", color="dodgerblue", s=20, zorder=3)
	# ~ ax.scatter(0.1, -2.0, marker="x", color="dodgerblue", s=20, zorder=3)
	# ~ ax.scatter(0.3, -4.25, marker="x", color="dodgerblue", s=20, zorder=3)
	# ~ ax.scatter(0.37, -4.1, marker="x", color="dodgerblue", s=20, zorder=3)
	# ~ ax.scatter(0.3, -3.8, marker="x", color="dodgerblue", s=20, zorder=3)
	# ~ ax.text(0.31, 2.95, "1)", fontsize=12)
	# ~ ax.text(0.11, -2.15, "2)", fontsize=12)
	# ~ ax.text(0.31, -4.4, "3)", fontsize=12)
	# ~ ax.text(0.38, -3.9, "4)", fontsize=12)
	# ~ ax.text(0.31, -3.85, "5)", color="white", fontsize=12)
	
	# period diagram ###################################################
	
	# plot 2D injection period diagram
	fig, ax = plt.subplots()
	fig.set_size_inches(5.9, 5.0)
	fig.subplots_adjust(top=1.07, bottom=0.04, left=0.15, right=0.79)
	# ~ plt.suptitle(suptitle + "period")
	# ~ plt.title(title)
	period = plt.imshow(periodData, origin="lower", extent=extent, norm=mpl.colors.LogNorm(), cmap="magma", vmin=1e-1, vmax=1e1)
	chaos  = plt.imshow(chaoticRegion, origin="lower", extent=extent, cmap=mycolormap2)
	# ~ plt.plot([K_inj_line, K_inj_line], [-6.0, 6.0], color="black")
	plt.xlabel(xlabel, fontsize=18)
	plt.ylabel(ylabel, fontsize=18)
	plt.xticks(xticks, xticksNames, fontsize=18)
	plt.yticks(fontsize=18)
	ax.set_aspect(aspect)
	
	ax2 = fig.add_axes([0.82, 0.40, 0.03, 0.56])
	ax3 = fig.add_axes([0.82, 0.15, 0.03, 0.25])
	
	# divider = make_axes_locatable(ax)
	# cax = divider.append_axes("right", size="5%", pad=0.1)
	# cb = plt.colorbar(fraction=0.046, pad=0.04)
	
	cb_period = plt.colorbar(period, ax2)
	cb_period.set_label(r"period $T$ / ns", fontsize=18, labelpad=-5.0)
	cb_period.ax.tick_params(labelsize=18)
	
	cb_chaos = plt.colorbar(chaos, ax3)
	cb_chaos.set_label("complex", fontsize=18, labelpad=44.0)
	cb_chaos.ax.set_yticklabels([])
	cb_chaos.set_ticks([])
	
	# ~ ax.scatter(0.1, -1.2588, color="black", s=20)
	# ~ ax.scatter(0.3, 2.647, color="red", s=20, zorder=3)
	# ~ ax.scatter(0.3, 2.8, marker="x", color="dodgerblue", s=20, zorder=3)
	# ~ ax.scatter(0.1, -2.0, marker="x", color="dodgerblue", s=20, zorder=3)
	# ~ ax.scatter(0.3, -4.25, marker="x", color="dodgerblue", s=20, zorder=3)
	# ~ ax.scatter(0.37, -4.1, marker="x", color="dodgerblue", s=20, zorder=3)
	# ~ ax.scatter(0.3, -3.8, marker="x", color="dodgerblue", s=20, zorder=3)
	# ~ ax.text(0.31, 2.95, "1)", color="white", fontsize=12)
	# ~ ax.text(0.11, -2.15, "2)", color="white", fontsize=12)
	# ~ ax.text(0.31, -4.4, "3)", color="white", fontsize=12)
	# ~ ax.text(0.38, -3.9, "4)", fontsize=12)
	# ~ ax.text(0.31, -3.85, "5)", fontsize=12)
	'''
	
	# Lyapunov exponent diagram ########################################
	
	fig, ax = plt.subplots()
	fig.set_size_inches(5.9, 5.0)
	fig.subplots_adjust(top=1.07, bottom=0.04, left=0.15, right=0.79)
	# ~ plt.suptitle(suptitle + "Ljapunov exponent")
	# ~ plt.title(title)
	LyapunovPos = plt.imshow(LyapunovDataPos, origin="lower", extent=extent, cmap="Reds")
	LyapunovNeg = plt.imshow(LyapunovDataNeg, origin="lower", extent=extent, cmap="winter")
	plt.xlabel(xlabel, fontsize=18)
	plt.ylabel(ylabel, fontsize=18)
	plt.xticks(xticks, xticksNames, fontsize=18)
	plt.yticks(fontsize=18)
	ax.set_aspect(aspect)
	
	# ~ ax.scatter(0.1, -1.8850417, color="black", s=20)
	
	ax2 = fig.add_axes([0.82, 0.555, 0.03, 0.405])
	ax3 = fig.add_axes([0.82, 0.15, 0.03, 0.405])
	
	cb_LyapunovPos = plt.colorbar(LyapunovPos, ax2)
	cb_LyapunovPos.set_label(r"greatest Lyapunov exponent $\Lambda_{max}$", y=-0.1, labelpad=15.0, fontsize=18)
	cb_LyapunovPos.ax.set_yticklabels([0])
	cb_LyapunovPos.ax.tick_params(size=0.0, labelsize=18)
	
	cb_LyapunovNeg = plt.colorbar(LyapunovNeg, ax3)
	cb_LyapunovNeg.ax.set_yticklabels([])
	cb_LyapunovNeg.ax.tick_params(size=0.0, labelsize=18)
	
	plt.show()



# plots 2D injection diagrams for MDRE model: network
def plotInjectionNetwork ():
	
	# simulation parameters
	Ttrans = 1000
	Teval  = 50
	beta   = 0
	K_fb   = 0.05
	phi    = 0
	tau    = 0.8
	
	# parameters for 257 x 385 map
	KinjStart   = 0.0
	KinjEnd     = 0.5
	KinjStep    = 0.001953125
	DnuInjStart = -6.0
	DnuInjEnd   = 6.0
	DnuInjStep  = 0.03125
	
	# initialize data and state matrix
	DnuInjSteps = int((DnuInjEnd - DnuInjStart) / DnuInjStep) + 1
	KinjSteps   = int((KinjEnd - KinjStart) / KinjStep) + 1
	
	cwIntensityData0 = np.zeros((DnuInjSteps, KinjSteps))
	maximaData0      = np.zeros((DnuInjSteps, KinjSteps))
	periodData0      = np.zeros((DnuInjSteps, KinjSteps))
	cwIntensityData1 = np.zeros((DnuInjSteps, KinjSteps))
	maximaData1      = np.zeros((DnuInjSteps, KinjSteps))
	periodData1      = np.zeros((DnuInjSteps, KinjSteps))
	cwIntensityData2 = np.zeros((DnuInjSteps, KinjSteps))
	maximaData2      = np.zeros((DnuInjSteps, KinjSteps))
	periodData2      = np.zeros((DnuInjSteps, KinjSteps))
	cwIntensityData3 = np.zeros((DnuInjSteps, KinjSteps))
	maximaData3      = np.zeros((DnuInjSteps, KinjSteps))
	periodData3      = np.zeros((DnuInjSteps, KinjSteps))
	deviation01_data = np.zeros((DnuInjSteps, KinjSteps))
	deviation02_data = np.zeros((DnuInjSteps, KinjSteps))
	deviation03_data = np.zeros((DnuInjSteps, KinjSteps))
	deviation12_data = np.zeros((DnuInjSteps, KinjSteps))
	deviation13_data = np.zeros((DnuInjSteps, KinjSteps))
	deviation23_data = np.zeros((DnuInjSteps, KinjSteps))
	states_data      = np.zeros((DnuInjSteps, KinjSteps))
	
	
	# read data and fill data matrix
	
	i = 0
	
	for K_inj in np.arange(KinjStart, KinjEnd + 1e-4, KinjStep):
		
		cwIntensity0 = []
		maxima0      = []
		period0      = []
		cwIntensity1 = []
		maxima1      = []
		period1      = []
		cwIntensity2 = []
		maxima2      = []
		period2      = []
		cwIntensity3 = []
		maxima3      = []
		period3      = []
		deviation01  = []
		deviation02  = []
		deviation03  = []
		deviation12  = []
		deviation13  = []
		deviation23  = []
		
		with open("four_node_all_to_all/tau=8e-2/K_fb=0,05_phi=0,5pi_tau=8e-2/injection_line_scan_mdre_K_inj={:f}".format(K_inj) + ".txt", "r") as file:
			lines = file.readlines()
		
		for line in lines:
			cwIntensity0.append(float((line.split('	')[2])))
			maxima0.append(float((line.split('	')[3])))
			period0.append(float((line.split('	')[5])))
			cwIntensity1.append(float((line.split('	')[6])))
			maxima1.append(float((line.split('	')[7])))
			period1.append(float((line.split('	')[9])))
			cwIntensity2.append(float((line.split('	')[10])))
			maxima2.append(float((line.split('	')[11])))
			period2.append(float((line.split('	')[13])))
			cwIntensity3.append(float((line.split('	')[14])))
			maxima3.append(float((line.split('	')[15])))
			period3.append(float((line.split('	')[17])))
			deviation01.append(float((line.split('	')[18])))
			deviation02.append(float((line.split('	')[19])))
			deviation03.append(float((line.split('	')[20])))
			deviation12.append(float((line.split('	')[21])))
			deviation13.append(float((line.split('	')[22])))
			deviation23.append(float((line.split('	')[23])))
		
		cwIntensity0 = np.array(cwIntensity0)
		maxima0      = np.array(maxima0)
		period0      = np.array(period0)
		cwIntensity1 = np.array(cwIntensity1)
		maxima1      = np.array(maxima1)
		period1      = np.array(period1)
		cwIntensity2 = np.array(cwIntensity2)
		maxima2      = np.array(maxima2)
		period2      = np.array(period2)
		cwIntensity3 = np.array(cwIntensity3)
		maxima3      = np.array(maxima3)
		period3      = np.array(period3)
		deviation01  = np.array(deviation01)
		deviation02  = np.array(deviation02)
		deviation03  = np.array(deviation03)
		deviation12  = np.array(deviation12)
		deviation13  = np.array(deviation13)
		deviation23  = np.array(deviation23)
		
		cwIntensityData0[:, i] = cwIntensity0
		maximaData0[:, i]      = maxima0
		periodData0[:, i]      = period0
		cwIntensityData1[:, i] = cwIntensity1
		maximaData1[:, i]      = maxima1
		periodData1[:, i]      = period1
		cwIntensityData2[:, i] = cwIntensity2
		maximaData2[:, i]      = maxima2
		periodData2[:, i]      = period2
		cwIntensityData3[:, i] = cwIntensity3
		maximaData3[:, i]      = maxima3
		periodData3[:, i]      = period3
		deviation01_data[:, i] = deviation01
		deviation02_data[:, i] = deviation02
		deviation03_data[:, i] = deviation03
		deviation12_data[:, i] = deviation12
		deviation13_data[:, i] = deviation13
		deviation23_data[:, i] = deviation23
		
		i += 1
	
	
	# general title and labels
	suptitle = "MDRE - four node all-to-all: injection - "
	title    = r"$\beta=$" + str(beta) + r", $K_{fb}=$" + str(K_fb) + r", $\varphi=\frac{3\pi}{2}$" + r", $\tau=$" + str(tau) + r"$\,$ns, $T_{trans}=$" + str(Ttrans) + r"$\,$ns, $T_{eval}=$" + str(Teval) + r"$\,$ns"
	xlabel      = r"injection strength $K_{inj}$"
	ylabel      = r"master-slave detuning $\Delta\nu_{inj}$ / GHz"		#$\longrightarrow$
	xticks      = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
	xticksNames = ["0", "0.1", "0.2", "0.3", "0.4", "0.5"]
	
	
	# bound maxima
	maximaData0[np.where(maximaData0 > 4)] = 5
	maximaData1[np.where(maximaData1 > 4)] = 5
	maximaData2[np.where(maximaData2 > 4)] = 5
	maximaData3[np.where(maximaData3 > 4)] = 5
	
	
	# define custom colors for ..
	
	# .. maxima plot
	mylightyellow = "#fff9bd"
	mydarkyellow  = "#fecb67"
	myorange      = "#fd8a3b"
	myred         = "#e41c1d"
	mykarmin      = "#be0126"
	
	# .. chaos regions in period plot
	mygray        = "#bbbbbb"
	
	# .. network states plot
	lexicon       = "#F2F2F2"
	coral         = "#F2B96E"
	melon         = "#E48251"
	redCapital    = "#9E3329"
	seaSight      = "#27839C"
	deepArctic    = "#264A4A"
	
	
	# define custom color maps, bounds and ticks
	
	mycolormap     = mpl.colors.ListedColormap([mylightyellow, mydarkyellow, myorange, myred, mykarmin])
	colorBarTicks  = np.arange(5)
	colorBarBounds = np.arange(6) - 0.5
	
	mycolormap2    = mpl.colors.ListedColormap([mygray])
	
	mycolormap3    = mpl.colors.ListedColormap([lexicon, coral, melon, redCapital, seaSight, deepArctic])
	colorBarTicks3 = np.arange(6)
	colorBarBound3 = np.arange(7) - 0.5
	
	extent = np.array([KinjStart, KinjEnd, DnuInjStart, DnuInjEnd])
	aspect = (KinjEnd - KinjStart) / (DnuInjEnd - DnuInjStart)
	
	
	# big magic is happening here
	mappable  = mpl.cm.ScalarMappable(cmap=mycolormap)
	mappable.set_array(np.arange(5))
	mappable3 = mpl.cm.ScalarMappable(cmap=mycolormap3)
	mappable3.set_array(np.arange(6))
	
	
	# masking the arrays
	
	maximaDataBool   = maximaData0 == 0
	maximaData0      = np.ma.masked_where(maximaDataBool > 0, maximaData0)
	cwIntensityData0 = np.ma.masked_where(maximaDataBool == 0, cwIntensityData0)
	
	maximaDataBool   = maximaData1 == 0
	maximaData1      = np.ma.masked_where(maximaDataBool > 0, maximaData1)
	cwIntensityData1 = np.ma.masked_where(maximaDataBool == 0, cwIntensityData1)
	
	maximaDataBool   = maximaData2 == 0
	maximaData2      = np.ma.masked_where(maximaDataBool > 0, maximaData2)
	cwIntensityData2 = np.ma.masked_where(maximaDataBool == 0, cwIntensityData2)
	
	maximaDataBool   = maximaData3 == 0
	maximaData3      = np.ma.masked_where(maximaDataBool > 0, maximaData3)
	cwIntensityData3 = np.ma.masked_where(maximaDataBool == 0, cwIntensityData3)
	
	maximaDataBool   = maximaData0 == 5
	chaoticRegion0   = np.ones(periodData0.shape)
	periodData0      = np.ma.masked_where(maximaDataBool == True, periodData0)
	chaoticRegion0   = np.ma.masked_where(maximaDataBool == False, chaoticRegion0)
	
	maximaDataBool   = maximaData1 == 5
	chaoticRegion1   = np.ones(periodData1.shape)
	periodData1      = np.ma.masked_where(maximaDataBool == True, periodData1)
	chaoticRegion1   = np.ma.masked_where(maximaDataBool == False, chaoticRegion1)
	
	maximaDataBool   = maximaData2 == 5
	chaoticRegion2   = np.ones(periodData2.shape)
	periodData2      = np.ma.masked_where(maximaDataBool == True, periodData2)
	chaoticRegion2   = np.ma.masked_where(maximaDataBool == False, chaoticRegion2)
	
	maximaDataBool   = maximaData3 == 5
	chaoticRegion3   = np.ones(periodData3.shape)
	periodData3      = np.ma.masked_where(maximaDataBool == True, periodData3)
	chaoticRegion3   = np.ma.masked_where(maximaDataBool == False, chaoticRegion3)
	
	
	# network states
	
	synch_threshold = 1e-6
	splay_threshold = 1e-3
	
	# synchronization states
	d01_s = deviation01_data <= synch_threshold
	d02_s = deviation02_data <= synch_threshold
	d03_s = deviation03_data <= synch_threshold
	d12_s = deviation12_data <= synch_threshold
	d13_s = deviation13_data <= synch_threshold
	d23_s = deviation23_data <= synch_threshold
	
	# asynchron states
	notd01_s = np.logical_not(d01_s)
	notd02_s = np.logical_not(d02_s)
	notd03_s = np.logical_not(d03_s)
	notd12_s = np.logical_not(d12_s)
	notd13_s = np.logical_not(d13_s)
	notd23_s = np.logical_not(d23_s)
	
	# differences
	d0102_d = np.absolute(deviation01_data - deviation02_data) < splay_threshold
	d0103_d = np.absolute(deviation01_data - deviation03_data) < splay_threshold
	d0112_d = np.absolute(deviation01_data - deviation12_data) < splay_threshold
	d0113_d = np.absolute(deviation01_data - deviation13_data) < splay_threshold
	d0123_d = np.absolute(deviation01_data - deviation23_data) < splay_threshold
	d0203_d = np.absolute(deviation02_data - deviation03_data) < splay_threshold
	d0212_d = np.absolute(deviation02_data - deviation12_data) < splay_threshold
	d0213_d = np.absolute(deviation02_data - deviation13_data) < splay_threshold
	d0223_d = np.absolute(deviation02_data - deviation23_data) < splay_threshold
	d0312_d = np.absolute(deviation03_data - deviation12_data) < splay_threshold
	d0313_d = np.absolute(deviation03_data - deviation13_data) < splay_threshold
	d0323_d = np.absolute(deviation03_data - deviation23_data) < splay_threshold
	d1213_d = np.absolute(deviation12_data - deviation13_data) < splay_threshold
	d1223_d = np.absolute(deviation12_data - deviation23_data) < splay_threshold
	d1323_d = np.absolute(deviation13_data - deviation23_data) < splay_threshold
	
	notd0102_d = np.logical_not(d0102_d)
	notd0103_d = np.logical_not(d0103_d)
	notd0112_d = np.logical_not(d0112_d)
	notd0113_d = np.logical_not(d0113_d)
	notd0123_d = np.logical_not(d0123_d)
	notd0203_d = np.logical_not(d0203_d)
	notd0212_d = np.logical_not(d0212_d)
	notd0213_d = np.logical_not(d0213_d)
	notd0223_d = np.logical_not(d0223_d)
	notd0312_d = np.logical_not(d0312_d)
	notd0313_d = np.logical_not(d0313_d)
	notd0323_d = np.logical_not(d0323_d)
	notd1213_d = np.logical_not(d1213_d)
	notd1223_d = np.logical_not(d1223_d)
	notd1323_d = np.logical_not(d1323_d)
	
	
	# cluster 3-1 states
	d010212_c31 = np.logical_and.reduce([d01_s, d02_s, d12_s, notd03_s, notd13_s, notd23_s])
	d121323_c31 = np.logical_and.reduce([d12_s, d13_s, d23_s, notd01_s, notd02_s, notd03_s])
	d020323_c31 = np.logical_and.reduce([d02_s, d03_s, d23_s, notd01_s, notd12_s, notd13_s])
	d010313_c31 = np.logical_and.reduce([d01_s, d03_s, d13_s, notd02_s, notd12_s, notd23_s])
	
	# cluster 2-2 states
	d0123_c22 = np.logical_and.reduce([d01_s, d23_s, notd02_s, notd03_s, notd12_s, notd13_s])
	d0213_c22 = np.logical_and.reduce([d02_s, d13_s, notd01_s, notd03_s, notd12_s, notd23_s])
	d0312_c22 = np.logical_and.reduce([d03_s, d12_s, notd01_s, notd02_s, notd13_s, notd23_s])
	
	# chimera states
	d01_c = np.logical_and.reduce([d01_s, notd02_s, notd03_s, notd12_s, notd13_s, notd23_s])
	d02_c = np.logical_and.reduce([d02_s, notd01_s, notd03_s, notd12_s, notd13_s, notd23_s])
	d03_c = np.logical_and.reduce([d03_s, notd01_s, notd02_s, notd12_s, notd13_s, notd23_s])
	d12_c = np.logical_and.reduce([d12_s, notd01_s, notd02_s, notd03_s, notd13_s, notd23_s])
	d13_c = np.logical_and.reduce([d13_s, notd01_s, notd02_s, notd03_s, notd12_s, notd23_s])
	d23_c = np.logical_and.reduce([d23_s, notd01_s, notd02_s, notd03_s, notd12_s, notd13_s])
	
	# splay states
	asynchron = np.logical_and.reduce([notd01_s, notd02_s, notd03_s, notd12_s, notd13_s, notd23_s])
	d0123_sp  = np.logical_and.reduce([asynchron, d0123_d, d0203_d, d0212_d, d0213_d, d0312_d, d0313_d, d1213_d, notd0102_d, notd0103_d, notd0112_d, notd0113_d, notd0223_d, notd0323_d, notd1223_d, notd1323_d])
	d0213_sp  = np.logical_and.reduce([asynchron, d0103_d, d0112_d, d0123_d, d0213_d, d0312_d, d0323_d, d1223_d, notd0102_d, notd0113_d, notd0203_d, notd0212_d, notd0223_d, notd0313_d, notd1213_d, notd1323_d])
	d0312_sp  = np.logical_and.reduce([asynchron, d0102_d, d0113_d, d0123_d, d0213_d, d0223_d, d0312_d, d1323_d, notd0103_d, notd0112_d, notd0203_d, notd0212_d, notd0313_d, notd0323_d, notd1213_d, notd1223_d])
	
	synchronized   = np.logical_and.reduce([d01_s, d02_s, d03_s, d12_s, d13_s, d23_s])
	cluster31      = np.logical_or.reduce([d010212_c31, d121323_c31, d020323_c31, d010313_c31])
	cluster22      = np.logical_or.reduce([d0123_c22, d0213_c22, d0312_c22])
	chimera        = np.logical_or.reduce([d01_c, d02_c, d03_c, d12_c, d13_c, d23_c])
	splay          = np.logical_or.reduce([d0123_sp, d0213_sp, d0312_sp])
	desynchronized = np.logical_not(np.logical_or.reduce([synchronized, cluster31, cluster22, chimera, splay]))
	
	states_data[np.where(synchronized == 1)]   = 1
	states_data[np.where(cluster31 == 1)]      = 2
	states_data[np.where(cluster22 == 1)]      = 3
	states_data[np.where(chimera == 1)]        = 4
	states_data[np.where(splay == 1)]          = 5
	states_data[np.where(desynchronized == 1)] = 6
	
	
	# maxima diagram ###################################################
	
	fig, ax = plt.subplots(1, 1)
	fig.set_size_inches(5.9, 5.0)
	fig.subplots_adjust(top=1.07, bottom=0.04, left=0.15, right=0.79)
	# ~ plt.suptitle(suptitle + "unique maxima", fontsize=14)
	# ~ plt.title(title)
	
	maxima0 = ax.imshow(maximaData0, origin="lower", extent=extent, cmap=mycolormap)
	cwInt0  = ax.imshow(cwIntensityData0, origin="lower", extent=extent, cmap="Blues_r")
	
	plt.xlabel(xlabel, fontsize=18)
	plt.ylabel(ylabel, fontsize=18)
	plt.xticks(xticks, xticksNames, fontsize=18)
	plt.yticks(fontsize=18)
	ax.set_aspect(aspect)
	
	# ~ ax.scatter(0.1, -1.8850417, color="black", s=20)
	
	ax2 = fig.add_axes([0.82, 0.575, 0.03, 0.385])
	ax3 = fig.add_axes([0.82, 0.15, 0.03, 0.385])
	
	cb_int = fig.colorbar(cwInt0, ax2)
	cb_int.set_label(r"cw intensity", fontsize=18, labelpad=2.0) # 2.0 for phi=0, 1.5pi | 19.0 for phi=0.5pi, pi | 
	cb_int.ax.tick_params(labelsize=18)
	
	cb_max0 = fig.colorbar(mappable, ax3, boundaries=colorBarBounds, ticks=colorBarTicks)
	cb_max0.ax.set_yticklabels(["1", "2", "3", "4", ">4"], fontsize=18)
	cb_max0.set_label("unique maxima", fontsize=18)
	
	
	# period diagram ###################################################
	
	fig, ax = plt.subplots(1, 1)
	fig.set_size_inches(5.9, 5.0)
	fig.subplots_adjust(top=1.07, bottom=0.04, left=0.15, right=0.79)
	# ~ plt.suptitle(suptitle + "period", fontsize=14)
	# ~ plt.title(title)
	
	period = ax.imshow(periodData0, origin="lower", extent=extent, norm=mpl.colors.LogNorm(), cmap="magma", vmin=1e-1, vmax=1e1)
	chaos  = plt.imshow(chaoticRegion0, origin="lower", extent=extent, cmap=mycolormap2)
	
	plt.xlabel(xlabel, fontsize=18)
	plt.ylabel(ylabel, fontsize=18)
	plt.xticks(xticks, xticksNames, fontsize=18)
	plt.yticks(fontsize=18)
	ax.set_aspect(aspect)
	
	# ~ ax.scatter(0.1, -1.8850417, color="black", s=20)
	
	ax2 = fig.add_axes([0.82, 0.40, 0.03, 0.56])
	ax3 = fig.add_axes([0.82, 0.15, 0.03, 0.25])
	
	cb_period = plt.colorbar(period, ax2)
	cb_period.set_label(r"period $T$ / ns", fontsize=18, labelpad=-2.5)
	cb_period.ax.tick_params(labelsize=18)
	
	cb_chaos = plt.colorbar(chaos, ax3)
	cb_chaos.set_label("complex", fontsize=18, labelpad=47.0)
	cb_chaos.ax.set_yticklabels([])
	cb_chaos.set_ticks([])
	
	
	# network states diagram ###########################################
	
	fig, ax = plt.subplots(1, 1)
	fig.set_size_inches(5.9, 5.0)
	fig.subplots_adjust(top=1.07, bottom=0.04, left=0.15, right=0.79)
	# ~ plt.suptitle(suptitle + "network states", fontsize=14)
	# ~ plt.title(title)
	
	states = ax.imshow(states_data, origin="lower", extent=extent, cmap=mycolormap3)
	
	# ~ ax.scatter(0.1, -1.8850417, color="white", s=20)
	
	plt.xlabel(xlabel, fontsize=18)
	plt.ylabel(ylabel, fontsize=18)
	plt.xticks(xticks, xticksNames, fontsize=18)
	plt.yticks(fontsize=18)
	ax.set_aspect(aspect)
	
	# ~ ax.scatter(0.1, -0.8, color="black", s=10)
	# ~ ax.scatter(0.2,  4.0, color="black", s=10)
	# ~ ax.scatter(0.16, 1.9, color="black", s=10)
	# ~ ax.scatter(0.05, 0.75, color="black", s=10)
	# ~ ax.scatter(0.12, 2.0, color="white", s=10)
	# ~ ax.scatter(0.19, 2.4, color="black", s=10)
	# ~ ax.scatter(0.02, -1.0, color="black", s=10)
	# ~ ax.scatter(0.05, -3.0, color="black", s=10)
	# ~ ax.scatter(0.11, -3.0, color="white", s=10)
	# ~ ax.scatter(0.15, -3.4, color="white", s=10)
	# ~ ax.scatter(0.27, -3.67, color="black", s=10)
	# ~ ax.scatter(0.35, -5.5, color="black", s=10)
	
	ax2 = fig.add_axes([0.82, 0.15, 0.03, 0.81])
	
	cb = fig.colorbar(mappable3, ax2, boundaries=colorBarBound3, ticks=colorBarTicks3, pad=0.02)
	cb.ax.set_yticklabels(["synchr", "clust 13", "clust 22", "chim", "splay", "desyn"], fontsize=14)
	# ~ cb.ax.set_yticklabels(["synchronized", "cluster 3-1", "cluster 2-2", "chimera", "splay", "desynchronized"])
	# cb.set_label("network states")
	
	plt.show()


# ~ plotInjectionSingleLaser()
plotInjectionNetwork()


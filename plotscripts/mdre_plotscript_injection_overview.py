import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


def readData (path):
	
	# parameters for 257 x 257 map
	K_injStart  = 0.0
	K_injEnd    = 0.5
	K_injStep   = 0.001953125
	DnuInjStart = -6.0
	DnuInjEnd   = 6.0
	DnuInjStep  = 0.03125
	
	# initialize data and matrices
	DnuInjSteps = int((DnuInjEnd - DnuInjStart) / DnuInjStep) + 1
	K_injSteps  = int((K_injEnd - K_injStart) / K_injStep) + 1
	
	cwIntensityData = np.zeros((DnuInjSteps, K_injSteps))
	maximaData      = np.zeros((DnuInjSteps, K_injSteps))
	
	i = 0
	
	for K_inj in np.arange(K_injStart, K_injEnd + 1e-4, K_injStep):
		
		cwIntensity = []
		maxima      = []
		
		with open(path + "injection_line_scan_mdre_K_inj={:f}".format(K_inj) + ".txt", "r") as file:
			lines = file.readlines()
		
		for line in lines:
			cwIntensity.append(float((line.split('	')[2])))
			maxima.append(int((line.split('	')[3])))
		
		cwIntensity = np.array(cwIntensity)
		maxima      = np.array(maxima)
		
		cwIntensityData[:, i] = cwIntensity
		maximaData[:, i]      = maxima
		
		i += 1
	
	# bound maxima
	maximaData[np.where(maximaData > 4)] = 5
	
	# masking the arrays
	maximaDataBool  = maximaData == 0
	maximaData      = np.ma.masked_where(maximaDataBool > 0, maximaData)
	cwIntensityData = np.ma.masked_where(maximaDataBool == 0, cwIntensityData)
	
	return (cwIntensityData, maximaData)


# plots 2D injection diagram overview for MDRE model: single laser
def plotInjectionSingleLaserOverview ():
	
	# simulation parameters
	Ttrans = 1000
	Teval  = 50
	K_fb   = 0.05
	phi    = 0.0
	tau    = 0.8
	
	# parameters for 257 x 385 map
	KinjStart   = 0.0
	KinjEnd     = 0.5
	KinjStep    = 0.001953125
	DnuInjStart = -6.0
	DnuInjEnd   = 6.0
	DnuInjStep  = 0.03125
	
	K_inj_line = 0.3
	
	# initialize data and state matrix
	DnuInjSteps = int((DnuInjEnd - DnuInjStart) / DnuInjStep) + 1
	KinjSteps   = int((KinjEnd - KinjStart) / KinjStep) + 1
	
	cwIntensityData = np.zeros((DnuInjSteps, KinjSteps))
	maximaData      = np.zeros((DnuInjSteps, KinjSteps))
	
	
	# pathlists for tau = 0.08 ns
	pathListSingleStrength_phi0_tau008 = ["single_Laser/tau=8e-2/K_fb=0,05_phi=0_tau=8e-2/", 
						"single_Laser/tau=8e-2/K_fb=0,15_phi=0_tau=8e-2/", 
						"single_Laser/tau=8e-2/K_fb=0,25_phi=0_tau=8e-2/", 
						"single_Laser/tau=8e-2/K_fb=0,35_phi=0_tau=8e-2/"]
	
	# ~ pathListSingleStrength_phipi_tau008 = ["single_Laser/tau=8e-2/K_fb=0,05_phi=pi_tau=8e-2/", 
						# ~ "single_Laser/tau=8e-2/K_fb=0,15_phi=pi_tau=8e-2/", 
						# ~ "single_Laser/tau=8e-2/K_fb=0,25_phi=pi_tau=8e-2/", 
						# ~ ""]
	
	pathListSinglePhase_tau008 = ["single_Laser/tau=8e-2/K_fb=0,05_phi=0_tau=8e-2/", 
				"single_Laser/tau=8e-2/K_fb=0,05_phi=0,5pi_tau=8e-2/", 
				"single_Laser/tau=8e-2/K_fb=0,05_phi=pi_tau=8e-2/", 
				"single_Laser/tau=8e-2/K_fb=0,05_phi=1,5pi_tau=8e-2/"]
	
	# pathlists for tau = 0.3 ns
	pathListSingleStrength_phi0_tau03 = ["single_Laser/tau=3e-1/K_fb=0,05_phi=0_tau=3e-1/", 
				"single_Laser/tau=3e-1/K_fb=0,15_phi=0_tau=3e-1/", 
				"single_Laser/tau=3e-1/K_fb=0,25_phi=0_tau=3e-1/", 
				"single_Laser/tau=3e-1/K_fb=0,35_phi=0_tau=3e-1/"]
	
	pathListSinglePhase_tau03 = ["single_Laser/tau=3e-1/K_fb=0,05_phi=0_tau=3e-1/", 
				"single_Laser/tau=3e-1/K_fb=0,05_phi=0,5pi_tau=3e-1/", 
				"single_Laser/tau=3e-1/K_fb=0,05_phi=pi_tau=3e-1/", 
				"single_Laser/tau=3e-1/K_fb=0,05_phi=1,5pi_tau=3e-1/"]
	
	# pathlists for tau = 0.8 ns
	pathListSingleStrength_phi0_tau08 = ["single_Laser/tau=8e-1/K_fb=0,05_phi=0_tau=8e-1/", 
				"single_Laser/tau=8e-1/K_fb=0,15_phi=0_tau=8e-1/", 
				"single_Laser/tau=8e-1/K_fb=0,25_phi=0_tau=8e-1/", 
				"single_Laser/tau=8e-1/K_fb=0,35_phi=0_tau=8e-1/"]
	
	pathListSinglePhase_tau08 = ["single_Laser/tau=8e-1/K_fb=0,05_phi=0_tau=8e-1/", 
				"single_Laser/tau=8e-1/K_fb=0,05_phi=0,5pi_tau=8e-1/", 
				"single_Laser/tau=8e-1/K_fb=0,05_phi=pi_tau=8e-1/", 
				"single_Laser/tau=8e-1/K_fb=0,05_phi=1,5pi_tau=8e-1/"]
	
	
	# read data and fill data matrix
	
	cwList  = []
	maxList = []
	
	i = 0
	
	# change paths here:
	for path in pathListSinglePhase_tau008:
		
		cwIntensityData, maximaData = readData(path)
		
		cwList.append(cwIntensityData)
		maxList.append(maximaData)
	
	
	# general title and labels
	xlabel      = r"injection strength $K_{inj}$"
	ylabel      = r"master-slave detuning $\Delta\nu_{inj}$ / GHz"		#$\longrightarrow$
	xticks      = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
	xticksNames = ["0", "0.1", "0.2", "0.3", "0.4", "0.5"]
	yticks      = [-6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0]
	yticksNames = ["-6", "-4", "-2", "0", "2", "4", "6"]
	
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
	
	
	# overview diagram #################################################
	
	fig, axes = plt.subplots(2, 2, sharex=True, sharey=True)
	ax        = fig.add_subplot(111, frameon=False)
	fig.set_size_inches(5.9, 4.9)
	fig.subplots_adjust(top=0.98, bottom=0.11, left=0.09, right=0.77)
	plt.rcParams.update({"font.size": 9})
	
	max00 = axes[0, 0].imshow(maxList[0], origin="lower", extent=extent, cmap=mycolormap)
	cw00  = axes[0, 0].imshow(cwList[0], origin="lower", extent=extent, cmap="Blues_r")
	axes[0, 0].set_aspect(aspect)
	
	max01 = axes[0, 1].imshow(maxList[1], origin="lower", extent=extent, cmap=mycolormap)
	cw01  = axes[0, 1].imshow(cwList[1], origin="lower", extent=extent, cmap="Blues_r")
	axes[0, 1].set_aspect(aspect)
	
	max10 = axes[1, 0].imshow(maxList[2], origin="lower", extent=extent, cmap=mycolormap)
	cw10  = axes[1, 0].imshow(cwList[2], origin="lower", extent=extent, cmap="Blues_r")
	axes[1, 0].set_aspect(aspect)
	
	max11 = axes[1, 1].imshow(maxList[3], origin="lower", extent=extent, cmap=mycolormap)
	cw11  = axes[1, 1].imshow(cwList[3], origin="lower", extent=extent, cmap="Blues_r")
	axes[1, 1].set_aspect(aspect)
	
	
	# texts for feedback coupling strength variation
	# ~ axes[0, 0].text(0.04, 0.91, r"$K_{fb}=0.05$", transform = axes[0, 0].transAxes)
	# ~ axes[0, 1].text(0.04, 0.91, r"$K_{fb}=0.15$", transform = axes[0, 1].transAxes)
	# ~ axes[1, 0].text(0.04, 0.91, r"$K_{fb}=0.25$", transform = axes[1, 0].transAxes)
	# ~ axes[1, 1].text(0.04, 0.91, r"$K_{fb}=0.35$", transform = axes[1, 1].transAxes)
	
	# texts for feedback coupling phase variation
	axes[0, 0].text(0.04, 0.91, r"$\varphi=0$", transform = axes[0, 0].transAxes)
	axes[0, 1].text(0.04, 0.91, r"$\varphi=\pi/2$", transform = axes[0, 1].transAxes)
	axes[1, 0].text(0.04, 0.91, r"$\varphi=\pi$", transform = axes[1, 0].transAxes)
	axes[1, 1].text(0.04, 0.91, r"$\varphi=3\pi/2$", transform = axes[1, 1].transAxes)
	
	plt.setp(ax.get_yticklabels(), visible=False)
	plt.setp(ax.get_xticklabels(), visible=False)
	ax.tick_params(axis="both", which="both", length=0.0)
	
	plt.setp(axes[0, 0], xticks=xticks, xticklabels=xticksNames, yticks=yticks)
	axes[0, 0].tick_params(labelsize=9)
	axes[1, 0].tick_params(labelsize=9)
	axes[1, 1].tick_params(labelsize=9)
	
	plt.xlabel(xlabel, fontsize=9.0, labelpad=24.0)
	plt.ylabel(ylabel, fontsize=9.0, labelpad=24.0)
	
	# ~ axes[1, 0].plot([K_inj_line, K_inj_line], [-6.0, 6.0], color="black", linewidth=0.75)
	axes[0, 1].scatter(0.1, -1.8850417, color="black", s=10)
	
	cbint_ax = fig.add_axes([0.805, 0.56, 0.03, 0.42])
	cbmax_ax = fig.add_axes([0.805, 0.11, 0.03, 0.42])
	
	cbint = fig.colorbar(mappable=cw00, cax=cbint_ax)
	cbint.set_label(r"cw intensity")
	
	cbmax = fig.colorbar(mappable=mappable, cax=cbmax_ax, boundaries=colorBarBounds, ticks=colorBarTicks)
	cbmax.ax.set_yticklabels(["1 maximum", "2 maxima", "3 maxima", "4 maxima", "chaos or\n>4 maxima"])
	
	plt.show()


plotInjectionSingleLaserOverview()


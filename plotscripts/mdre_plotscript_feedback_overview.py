import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


def readData (path):
	
	# parameters for 257 x 257 map
	K_fbStart = 0.0
	K_fbEnd   = 0.5
	K_fbStep  = 0.001953125
	phiStart  = 0.0
	phiEnd    = 2.0 * np.pi
	phiStep   = 2.0 * np.pi / 256.0
	
	# initialize data and matrices
	phiSteps       = int((phiEnd - phiStart) / phiStep) + 1
	K_fbSteps      = int((K_fbEnd - K_fbStart) / K_fbStep) + 1
	cwIntensityData = np.zeros((phiSteps, K_fbSteps))
	maximaData      = np.zeros((phiSteps, K_fbSteps))
	states_data     = np.zeros((phiSteps, K_fbSteps))
	deviation01_data = np.zeros((phiSteps, K_fbSteps))
	deviation02_data = np.zeros((phiSteps, K_fbSteps))
	deviation03_data = np.zeros((phiSteps, K_fbSteps))
	deviation12_data = np.zeros((phiSteps, K_fbSteps))
	deviation13_data = np.zeros((phiSteps, K_fbSteps))
	deviation23_data = np.zeros((phiSteps, K_fbSteps))
	
	i = 0
	
	for K_fb in np.arange(K_fbStart, K_fbEnd + 1e-4, K_fbStep):
		
		cwIntensity = []
		maxima      = []
		deviation01  = []
		deviation02  = []
		deviation03  = []
		deviation12  = []
		deviation13  = []
		deviation23  = []
		
		with open(path + "feedback_line_scan_mdre_K_fb={:f}".format(K_fb) + ".txt", "r") as file:
			lines = file.readlines()
		
		for line in lines:
			cwIntensity.append(float((line.split('	')[2])))
			maxima.append(float((line.split('	')[3])))
			deviation01.append(float((line.split('	')[18])))
			deviation02.append(float((line.split('	')[19])))
			deviation03.append(float((line.split('	')[20])))
			deviation12.append(float((line.split('	')[21])))
			deviation13.append(float((line.split('	')[22])))
			deviation23.append(float((line.split('	')[23])))
		
		cwIntensity = np.array(cwIntensity)
		maxima      = np.array(maxima)
		deviation01  = np.array(deviation01)
		deviation02  = np.array(deviation02)
		deviation03  = np.array(deviation03)
		deviation12  = np.array(deviation12)
		deviation13  = np.array(deviation13)
		deviation23  = np.array(deviation23)
		
		cwIntensityData[:, i] = cwIntensity
		maximaData[:, i]      = maxima
		deviation01_data[:, i] = deviation01
		deviation02_data[:, i] = deviation02
		deviation03_data[:, i] = deviation03
		deviation12_data[:, i] = deviation12
		deviation13_data[:, i] = deviation13
		deviation23_data[:, i] = deviation23
		
		i += 1
	
	# bound maxima
	maximaData[np.where(maximaData > 4)] = 5
	
	# masking the arrays
	maximaDataBool  = maximaData == 0
	maximaData      = np.ma.masked_where(maximaDataBool > 0, maximaData)
	cwIntensityData = np.ma.masked_where(maximaDataBool == 0, cwIntensityData)
	
	
	
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
	
	
	
	return (cwIntensityData, maximaData, states_data)



# plots 2D feedback coupling diagram for MDRE model: rotation of tear drop
def plotRotationOfTearDrop ():
	
	# parameters
	Ttrans = 1000
	Teval  = 50
	tau    = 0.08
	
	# parameters for 257 x 257 map
	K_fbStart = 0.0
	K_fbEnd   = 0.5
	K_fbStep  = 0.001953125
	phiStart  = 0.0
	phiEnd    = 2.0 * np.pi
	phiStep   = 2.0 * np.pi / 256.0
	
	# ~ K_fb_line = np.pi
	
	# initialize data matrices
	phiSteps        = int((phiEnd - phiStart) / phiStep) + 1
	K_fbSteps       = int((K_fbEnd - K_fbStart) / K_fbStep) + 1
	cwIntensityData = np.zeros((phiSteps, K_fbSteps))
	states_data     = np.zeros((phiSteps, K_fbSteps))
	
	pathListSingle = ["single_laser/tau=8e-2/", 
					"single_laser/tau=1e-1/", 
					"single_laser/tau=1,5e-1/", 
					"single_laser/tau=2e-1/", 
					"single_laser/tau=3e-1/scan/", 
					"single_laser/tau=4e-1/", 
					"single_laser/tau=5e-1/", 
					"single_laser/tau=6e-1/", 
					"single_laser/tau=8e-1/"]
	
	pathListNetwork = ["four_node_all_to_all/tau=8e-2/with_perturbation/", 
					"four_node_all_to_all/tau=1e-1/with_perturbation/", 
					"four_node_all_to_all/tau=1,5e-1/with_perturbation/", 
					"four_node_all_to_all/tau=2e-1/with_perturbation/", 
					"four_node_all_to_all/tau=3e-1/with_perturbation/", 
					"four_node_all_to_all/tau=4e-1/with_perturbation/", 
					"four_node_all_to_all/tau=5e-1/with_perturbation/", 
					"four_node_all_to_all/tau=6e-1/with_perturbation/", 
					"four_node_all_to_all/tau=8e-1/with_perturbation/"]
	
	# ~ pathListNetwork = ["four_node_all_to_all/tau=8e-2/without_perturbation/", 
					# ~ "four_node_all_to_all/tau=1e-1/without_perturbation/", 
					# ~ "four_node_all_to_all/tau=1,5e-1/without_perturbation/", 
					# ~ "four_node_all_to_all/tau=2e-1/without_perturbation/", 
					# ~ "four_node_all_to_all/tau=3e-1/without_perturbation/", 
					# ~ "four_node_all_to_all/tau=4e-1/without_perturbation/", 
					# ~ "four_node_all_to_all/tau=5e-1/without_perturbation/", 
					# ~ "four_node_all_to_all/tau=6e-1/without_perturbation/", 
					# ~ "four_node_all_to_all/tau=8e-1/without_perturbation/"]
	
	
	# read data and fill data matrix
	
	cwList  = []
	maxList = []
	stList  = []
	
	
	i = 0
	
	# change paths here:
	for path in pathListNetwork:
		
		cwIntensityData, maximaData, states_data = readData(path)
		
		cwList.append(cwIntensityData)
		maxList.append(maximaData)
		stList.append(states_data)

	
	# ~ maxs = np.zeros(9)
	# ~ mins = np.zeros(9)
	
	# ~ for i in range(len(cwList)):
		# ~ mins[i] = np.min(cwList[i])
	
	# ~ max = np.max(maxs)
	# ~ min = np.min(mins)
	
	
	# define custom colors
	mywhite       = "#ffffff"
	mylightyellow = "#fff9bd"
	mydarkyellow  = "#fecb67"
	myorange      = "#fd8a3b"
	myred         = "#e41c1d"
	mykarmin      = "#be0126"
	
	# .. network states plot
	lexicon       = "#F2F2F2"
	coral         = "#F2B96E"
	melon         = "#E48251"
	redCapital    = "#9E3329"
	seaSight      = "#27839C"
	deepArctic    = "#264A4A"
	
	# define custom color map, bounds and ticks
	
	mycolormap     = mpl.colors.ListedColormap([mylightyellow, mydarkyellow, myorange, myred, mykarmin])
	colorBarTicks  = np.arange(5)
	colorBarBounds = np.arange(6) - 0.5
	
	mycolormap3    = mpl.colors.ListedColormap([lexicon, coral, melon, redCapital, seaSight, deepArctic])
	colorBarTicks3 = np.arange(6)
	colorBarBound3 = np.arange(7) - 0.5
	
	extent      = [K_fbStart, K_fbEnd, phiStart, phiEnd]
	aspect      = K_fbEnd / phiEnd
	
	xticks      = [0.0, 0.25, 0.5]
	xticksNames = ["0", "0.25", "0.5"]
	yticks      = np.array([0.0, np.pi, 2.0 * np.pi])
	yticksNames = [0, r"$\pi$", r"$2\pi$"]
	
	# big magic is happening here
	mappable = mpl.cm.ScalarMappable(cmap=mycolormap)
	mappable.set_array(np.arange(5))
	mappable3 = mpl.cm.ScalarMappable(cmap=mycolormap3)
	mappable3.set_array(np.arange(6))
	
	
	fig, axes = plt.subplots(3, 3, sharex=True, sharey=True)
	fig.set_size_inches(5.9, 4.9)
	fig.subplots_adjust(top=0.98, bottom=0.11, left=0.09, right=0.77)
	plt.rcParams.update({"font.size": 9})
	
	# ~ max00 = axes[0, 0].imshow(maxList[0], origin="lower", extent=extent, cmap=mycolormap)
	# ~ cw00  = axes[0, 0].imshow(cwList[0], origin="lower", extent=extent, cmap="Blues_r")
	st00  = axes[0, 0].imshow(stList[0], origin="lower", extent=extent, cmap=mycolormap3)
	axes[0, 0].set_aspect(aspect)
	axes[0, 0].text(0.182, 5.44, r"$\tau=0.08\,$ns", bbox={"facecolor":"white", "edgecolor":"none", "alpha":0.5})
	
	# ~ max01 = axes[0, 1].imshow(maxList[1], origin="lower", extent=extent, cmap=mycolormap)
	# ~ cw01  = axes[0, 1].imshow(cwList[1], origin="lower", extent=extent, cmap="Blues_r")
	st01  = axes[0, 1].imshow(stList[1], origin="lower", extent=extent, cmap=mycolormap3)
	axes[0, 1].set_aspect(aspect)
	axes[0, 1].text(0.218, 5.44, r"$\tau=0.1\,$ns", bbox={"facecolor":"white", "edgecolor":"none", "alpha":0.5})
	
	# ~ max02 = axes[0, 2].imshow(maxList[2], origin="lower", extent=extent, cmap=mycolormap)
	# ~ cw01  = axes[0, 2].imshow(cwList[2], origin="lower", extent=extent, cmap="Blues_r")
	st01  = axes[0, 2].imshow(stList[2], origin="lower", extent=extent, cmap=mycolormap3)
	axes[0, 2].set_aspect(aspect)
	axes[0, 2].text(0.182, 5.44, r"$\tau=0.15\,$ns", bbox={"facecolor":"white", "edgecolor":"none", "alpha":0.5})
	
	# ~ max10 = axes[1, 0].imshow(maxList[3], origin="lower", extent=extent, cmap=mycolormap)
	# ~ cw10  = axes[1, 0].imshow(cwList[3], origin="lower", extent=extent, cmap="Blues_r")
	st10  = axes[1, 0].imshow(stList[3], origin="lower", extent=extent, cmap=mycolormap3)
	axes[1, 0].set_aspect(aspect)
	axes[1, 0].text(0.218, 5.44, r"$\tau=0.2\,$ns", bbox={"facecolor":"white", "edgecolor":"none", "alpha":0.5})
	
	# ~ max11 = axes[1, 1].imshow(maxList[4], origin="lower", extent=extent, cmap=mycolormap)
	# ~ cw11  = axes[1, 1].imshow(cwList[4], origin="lower", extent=extent, cmap="Blues_r")
	st11  = axes[1, 1].imshow(stList[4], origin="lower", extent=extent, cmap=mycolormap3)
	axes[1, 1].set_aspect(aspect)
	axes[1, 1].text(0.218, 5.44, r"$\tau=0.3\,$ns", bbox={"facecolor":"white", "edgecolor":"none", "alpha":0.5})
	
	# ~ max12 = axes[1, 2].imshow(maxList[5], origin="lower", extent=extent, cmap=mycolormap)
	# ~ cw12  = axes[1, 2].imshow(cwList[5], origin="lower", extent=extent, cmap="Blues_r")
	st12  = axes[1, 2].imshow(stList[5], origin="lower", extent=extent, cmap=mycolormap3)
	axes[1, 2].set_aspect(aspect)
	axes[1, 2].text(0.218, 5.44, r"$\tau=0.4\,$ns", bbox={"facecolor":"white", "edgecolor":"none", "alpha":0.5})
	
	# ~ max20 = axes[2, 0].imshow(maxList[6], origin="lower", extent=extent, cmap=mycolormap)
	# ~ cw20  = axes[2, 0].imshow(cwList[6], origin="lower", extent=extent, cmap="Blues_r")
	st20  = axes[2, 0].imshow(stList[6], origin="lower", extent=extent, cmap=mycolormap3)
	axes[2, 0].set_aspect(aspect)
	axes[2, 0].text(0.218, 5.44, r"$\tau=0.5\,$ns", bbox={"facecolor":"white", "edgecolor":"none", "alpha":0.5})
	
	# ~ max21 = axes[2, 1].imshow(maxList[7], origin="lower", extent=extent, cmap=mycolormap)
	# ~ cw21  = axes[2, 1].imshow(cwList[7], origin="lower", extent=extent, cmap="Blues_r")
	st21  = axes[2, 1].imshow(stList[7], origin="lower", extent=extent, cmap=mycolormap3)
	axes[2, 1].set_aspect(aspect)
	axes[2, 1].text(0.218, 5.44, r"$\tau=0.6\,$ns", bbox={"facecolor":"white", "edgecolor":"none", "alpha":0.5})
	
	# ~ max22 = axes[2, 2].imshow(maxList[8], origin="lower", extent=extent, cmap=mycolormap)
	# ~ cw22  = axes[2, 2].imshow(cwList[8], origin="lower", extent=extent, cmap="Blues_r")
	st22  = axes[2, 2].imshow(stList[8], origin="lower", extent=extent, cmap=mycolormap3)
	axes[2, 2].set_aspect(aspect)
	axes[2, 2].text(0.218, 5.44, r"$\tau=0.8\,$ns", bbox={"facecolor":"white", "edgecolor":"none", "alpha":0.5})
	
	plt.setp(axes[0, 0], xticks=xticks, xticklabels=xticksNames)
	plt.setp(axes[0, 0], yticks=yticks, yticklabels=yticksNames)
	
	axes[0, 0].tick_params(labelsize=9)
	axes[1, 0].tick_params(labelsize=9)
	axes[2, 0].tick_params(labelsize=9)
	axes[2, 1].tick_params(labelsize=9)
	axes[2, 2].tick_params(labelsize=9)
	
	axes[2, 1].set_xlabel(r"feedback coupling strength $K_{fb}$", fontsize=9)
	axes[1, 0].set_ylabel(r"feedback coupling phase $\varphi$", fontsize=9)
	
	# ~ cbint_ax = fig.add_axes([0.805, 0.56, 0.03, 0.42])
	# ~ cbmax_ax = fig.add_axes([0.805, 0.11, 0.03, 0.42])
	cbst_ax  = fig.add_axes([0.805, 0.11, 0.03, 0.87])
	
	# ~ cbint = fig.colorbar(mappable=cw00, cax=cbint_ax)
	# ~ cbint.set_label(r"cw intensity")
	
	# ~ cbmax = fig.colorbar(mappable=mappable, cax=cbmax_ax, boundaries=colorBarBounds, ticks=colorBarTicks)
	# ~ cbmax.ax.set_yticklabels(["1 maximum", "2 maxima", "3 maxima", "4 maxima", "chaos or\n>4 maxima"])
	
	cbst = fig.colorbar(mappable3, cax=cbst_ax, boundaries=colorBarBound3, ticks=colorBarTicks3, pad=0.05)
	cbst.ax.set_yticklabels(["synchronized", "cluster 1-3", "cluster 2-2", "chimera", "splay state", "desynchron."])
	
	plt.show()


plotRotationOfTearDrop()


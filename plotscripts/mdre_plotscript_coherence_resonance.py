import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


# plots coherence resonance curve for injection
def plotCoherenceResonanceInjection ():
	
	# parameters
	n_intervals = 1000
	Tint        = 350
	K_inj       = 0.1
	
	# parameters for 81 beta values
	betaExpStart = -3.5 + 4.0 * 0.03125
	betaExpEnd   = 0.0
	betaExpStep  = 0.03125
	
	# initialize data vectors
	betaExpSteps = int((betaExpEnd - betaExpStart) / betaExpStep) + 1
	betaData     = np.zeros(betaExpSteps)
	meanData     = np.zeros(betaExpSteps)
	varianceData = np.zeros(betaExpSteps)
	
	
	def parabolaFit (x, a, b, c):
		return (a * np.square(x) + b * x + c)
	
	
	# read data and fill data vectors
	
	i = 0
	
	for betaExp in np.arange(betaExpStart, betaExpEnd + 1e-4, betaExpStep):
		
		beta     = []
		mean     = []
		variance = []
		
		b = np.power(10.0, betaExp)
		
		with open("single_laser/injection/K_inj=0,1_n_intervals=1000_Tint=350_low=0,65_tdead=0,25_dt=1e-4_1e-15/coherence_resonance_beta={:f}".format(b) + ".txt", "r") as file:
			line = file.readline()
		
		beta.append(float((line.split('	')[0])))
		mean.append(float((line.split('	')[1])))
		variance.append(float((line.split('	')[2])))
		
		beta     = np.array(beta)
		mean     = np.array(mean)
		variance = np.array(variance)
		
		betaData[i]     = beta
		meanData[i]     = mean
		varianceData[i] = variance
		
		i += 1
	
	normalizedMean = meanData / np.max(meanData)
	normalizedStd  = np.sqrt(varianceData) / meanData
	
	
	# initial guess parameters for linear fit:
	a = 10.0
	b = 0.0
	c = 0.6
	initial_guess = [a, b, c]
	
	# fitting function parameters to data
	popt, pcov = curve_fit(parabolaFit, betaData, normalizedStd, p0=initial_guess)
	
	a = popt[0] #np.around(popt[0], 5)
	b = popt[1] #np.around(popt[1], 5)
	c = popt[2]
	
	a_err = pcov[0, 0] #np.around(np.sqrt(pcov[0, 0]), 6)
	b_err = pcov[1, 1] #np.around(np.sqrt(pcov[1, 1]), 6)
	c_err = pcov[2, 2] #np.around(np.sqrt(pcov[2, 2]), 6)
	
	
	plt.figure(figsize=(5.9, 4.6))
	plt.subplots_adjust(top=0.98, bottom=0.16, left=0.15, right=0.99)
	plt.rcParams.update({"font.size": 18})
	# ~ plt.suptitle("MDRE: coherence resonance", fontsize=14)
	# ~ plt.title(r"$K_{inj}=$" + str(K_inj) + r", $T_{int}=$" + str(Tint) +  r"$\,$ns, $n_{intervals}=$" + str(n_intervals))
	plt.semilogx(betaData, normalizedMean, color="red", marker="o", markersize=3, linewidth=0, label=r"$\mu_{T_{ISI}}}$")
	plt.semilogx(betaData, normalizedStd, color="darkblue", marker="o", markersize=3, linewidth=0, label=r"$\sigma_{T_{ISI}}$")
	# ~ plt.semilogx(betaData, parabolaFit(betaData, *popt))
	# ~ plt.xlim(2.75e-3, 1.1)
	# ~ plt.ylim(-0.025, np.max(np.sqrt(varianceData) / meanData) + 0.025)
	plt.xlabel(r"spontaneous emission coefficient $\beta$")
	plt.ylabel(r"cumulants of $T_{ISI,i}$")
	plt.axvline(x=1e-3, color="gray", linewidth=0.75)
	plt.axvline(x=1e-2, color="gray", linewidth=0.75)
	plt.axvline(x=1e-1, color="gray", linewidth=0.75)
	plt.legend(bbox_to_anchor=(1, 0.075), loc="lower right")
	plt.grid(color="lightgray")
	plt.show()


# plots coherence resonance curve for injection with delayed feedback
def plotCoherenceResonanceInjectionFeedback ():
	
	# parameters
	n_intervals = 1000
	Tint        = 1250
	K_inj       = 0.1
	
	# parameters for xx beta values
	betaExpStart = -2.5 + 2.0 * 0.03125
	betaExpEnd   = 0.0
	betaExpStep  = 0.03125
	
	# initialize data vectors
	betaExpSteps = int((betaExpEnd - betaExpStart) / betaExpStep) + 1
	betaData     = np.zeros(betaExpSteps)
	meanData     = np.zeros(betaExpSteps)
	varianceData = np.zeros(betaExpSteps)
	
	
	# read data and fill data vectors
	
	i = 0
	
	for betaExp in np.arange(betaExpStart, betaExpEnd + 1e-4, betaExpStep):
		
		beta     = []
		mean     = []
		variance = []
		
		b = np.power(10.0, betaExp)
		
		with open("single_laser/injection_feedback/K_fb=0,05_tau=8e-2_phi=0,5pi_K_inj=0,1_Tint=1250_low=0,70_tdead=0,25_dt=1e-4/coherence_resonance_beta={:f}".format(b) + ".txt", "r") as file:
			line = file.readline()
		
		beta.append(float((line.split('	')[0])))
		mean.append(float((line.split('	')[1])))
		variance.append(float((line.split('	')[2])))
		
		beta     = np.array(beta)
		mean     = np.array(mean)
		variance = np.array(variance)
		
		betaData[i]     = beta
		meanData[i]     = mean
		varianceData[i] = variance
		
		i += 1
	
	plt.figure(figsize=(5.9, 4.6))
	plt.subplots_adjust(top=0.98, bottom=0.16, left=0.15, right=0.99)
	plt.rcParams.update({"font.size": 18})
	# ~ plt.suptitle("MDRE: coherence resonance", fontsize=14)
	# ~ plt.title(r"$K_{inj}=$" + str(K_inj) + r", $T_{int}=$" + str(Tint) +  r"$\,$ns, $n_{intervals}=$" + str(n_intervals))
	plt.semilogx(betaData, meanData / np.max(meanData), color="red", marker="o", markersize=3, linewidth=0, label=r"$\mu_{T_{ISI}}}$")
	plt.semilogx(betaData, np.sqrt(varianceData) / meanData, color="darkblue", marker="o", markersize=3, linewidth=0, label=r"$\sigma_{T_{ISI}}$")
	# ~ plt.xlim(2.75e-3, 1.1)
	# ~ plt.ylim(-0.025, np.max(np.sqrt(varianceData) / meanData) + 0.025)
	plt.xlabel(r"spontaneous emission coefficient $\beta$")
	plt.ylabel(r"cumulants of $T_{ISI,i}$")
	plt.axvline(x=0.523299, color="gray", linewidth=0.75)
	plt.axvline(x=0.052330, color="gray", linewidth=0.75)
	plt.axvline(x=0.005233, color="gray", linewidth=0.75)
	
	plt.legend(bbox_to_anchor=(1, 0.075), loc="lower right")
	plt.grid(color="lightgray")
	plt.show()



# plots coherence resonance curve for four node all-to-all network
def plotCoherenceResonanceNetwork ():
	
	# parameters
	n_intervals = 1000
	Tint        = 1250
	K_inj       = 0.1
	tau         = 0.08
	phi         = np.pi / 2.0
	
	# parameters for xx beta values
	betaExpStart = -6.0#-5.5 + 15.0 * 0.03125
	betaExpEnd   = 0.0
	betaExpStep  = 0.03125
	
	# initialize data vectors
	betaExpSteps  = int((betaExpEnd - betaExpStart) / betaExpStep) + 1
	
	betaData      = np.zeros(betaExpSteps)
	meanData0     = np.zeros(betaExpSteps)
	meanData1     = np.zeros(betaExpSteps)
	meanData2     = np.zeros(betaExpSteps)
	meanData3     = np.zeros(betaExpSteps)
	varianceData0 = np.zeros(betaExpSteps)
	varianceData1 = np.zeros(betaExpSteps)
	varianceData2 = np.zeros(betaExpSteps)
	varianceData3 = np.zeros(betaExpSteps)
	
	
	# read data and fill data vectors
	
	i = 0
	
	for betaExp in np.arange(betaExpStart, betaExpEnd + 1e-4, betaExpStep):
		
		beta      = []
		mean0     = []
		mean1     = []
		mean2     = []
		mean3     = []
		variance0 = []
		variance1 = []
		variance2 = []
		variance3 = []
		
		b = np.power(10.0, betaExp)
		
		with open("network/coherence_resonance_beta={:f}".format(b) + ".txt", "r") as file:
			line = file.readline()
		
		beta.append(float((line.split('	')[0])))
		mean0.append(float((line.split('	')[1])))
		variance0.append(float((line.split('	')[2])))
		mean1.append(float((line.split('	')[3])))
		variance1.append(float((line.split('	')[4])))
		mean2.append(float((line.split('	')[5])))
		variance2.append(float((line.split('	')[6])))
		mean3.append(float((line.split('	')[7])))
		variance3.append(float((line.split('	')[8])))
		
		beta      = np.array(beta)
		mean0     = np.array(mean0)
		mean1     = np.array(mean1)
		mean2     = np.array(mean2)
		mean3     = np.array(mean3)
		variance0 = np.array(variance0)
		variance1 = np.array(variance1)
		variance2 = np.array(variance2)
		variance3 = np.array(variance3)
		
		betaData[i]      = beta
		meanData0[i]     = mean0
		varianceData0[i] = variance0
		meanData1[i]     = mean1
		varianceData1[i] = variance1
		meanData2[i]     = mean2
		varianceData2[i] = variance2
		meanData3[i]     = mean3
		varianceData3[i] = variance3
		
		i += 1
	
	plt.figure(figsize=(5.9, 4.6))
	plt.subplots_adjust(top=0.98, bottom=0.16, left=0.15, right=0.99)
	plt.rcParams.update({"font.size": 18})
	# ~ plt.suptitle("MDRE: coherence resonance", fontsize=14)
	# ~ plt.title(r"$K_{inj}=$" + str(K_inj) + r", $T_{int}=$" + str(Tint) +  r"$\,$ns, $n_{intervals}=$" + str(n_intervals))
	plt.semilogx(betaData, meanData0 / np.max(meanData0), color="red", marker="o", markersize=3, linewidth=0, label=r"$\mu_{T_{ISI,1}}}$")
	plt.semilogx(betaData, meanData1 / np.max(meanData1), color="red", marker=".", markersize=3, linewidth=0, label=r"$\mu_{T_{ISI,2}}}$")
	plt.semilogx(betaData, meanData2 / np.max(meanData2), color="red", marker="x", markersize=3, linewidth=0, label=r"$\mu_{T_{ISI,3}}}$")
	plt.semilogx(betaData, meanData3 / np.max(meanData3), color="red", marker="+", markersize=3, linewidth=0, label=r"$\mu_{T_{ISI,4}}}$")
	plt.semilogx(betaData, np.sqrt(varianceData0) / meanData0, color="darkblue", marker="o", markersize=3, linewidth=0, label=r"$\sigma_{T_{ISI,1}}$")
	plt.semilogx(betaData, np.sqrt(varianceData1) / meanData1, color="darkblue", marker=".", markersize=3, linewidth=0, label=r"$\sigma_{T_{ISI,2}}$")
	plt.semilogx(betaData, np.sqrt(varianceData2) / meanData2, color="darkblue", marker="x", markersize=3, linewidth=0, label=r"$\sigma_{T_{ISI,3}}$")
	plt.semilogx(betaData, np.sqrt(varianceData3) / meanData3, color="darkblue", marker="+", markersize=3, linewidth=0, label=r"$\sigma_{T_{ISI,4}}$")
	# ~ plt.xlim(2.75e-3, 1.1)
	# ~ plt.ylim(-0.025, np.max(np.sqrt(varianceData) / meanData) + 0.025)
	plt.xlabel(r"spontaneous emission coefficient $\beta$", fontsize=18)
	plt.ylabel(r"cumulants of $T_{ISI,i,j}$", fontsize=18)
	plt.xticks(fontsize=18)
	plt.yticks(fontsize=18)
	# ~ plt.legend(loc="lower right")
	plt.legend(bbox_to_anchor=(0.7, 0.235), loc="lower right", ncol=2, fontsize=18)
	plt.grid(color="lightgray")
	plt.show()



# plots histogram for arbitrary number of beta arguments
# for each beta a corresponding coherence resonance file  must exist 
def plotTISIHistogramInjection (*betaList):
	
	plt.figure(figsize=(5.9, 4.6))
	plt.subplots_adjust(top=0.99, bottom=0.17, left=0.17, right=0.96)
	plt.rcParams.update({"font.size": 18})
	
	binsLowLim = 1.0
	binsUpLim  = 1e3
	
	colorList = ["#00235E", "#618ABB", "#B9D1EA"]
	
	i = 0
	
	for b in betaList:
		# ~ TISI = np.loadtxt("single_laser/injection/coherence_resonance_beta={:f}".format(b) + ".txt", skiprows=2)
		TISI = np.loadtxt("single_laser/injection/K_inj=0,1_n_intervals=1000_Tint=350_low=0,65_tdead=0,25_dt=1e-4_1e-15/coherence_resonance_beta={:f}".format(b) + ".txt", skiprows=2)
		bins = np.logspace(np.log10(binsLowLim), np.log10(binsUpLim), 120)
		plt.hist(TISI, bins=bins, density=True, log=True, alpha=0.85, label=r"$\beta=$" + str(b), color=colorList[i])
		i += 1
	
	plt.xscale("log")
	plt.xlabel(r"interspike interval time $T_{ISI}$ / ns")
	plt.ylabel(r"frequency of $T_{ISI,i}$", labelpad=-1.0)
	plt.xlim(binsLowLim, binsUpLim)
	plt.legend(loc="upper right")
	plt.show()


# plots histogram for arbitrary number of beta arguments
# for each beta a corresponding coherence resonance file  must exist 
def plotTISIHistogramInjectionFeedback (*betaList):
	
	plt.figure(figsize=(5.9, 4.6))
	plt.subplots_adjust(top=0.99, bottom=0.17, left=0.17, right=0.96)
	plt.rcParams.update({"font.size": 18})
	
	binsLowLim = 0.1
	binsUpLim  = 1e3
	
	colorList = ["#00235E", "#618ABB", "#B9D1EA"]
	
	i = 0
	
	for b in betaList:
		TISI = np.loadtxt("single_laser/injection_feedback/K_fb=0,05_tau=8e-2_phi=0,5pi_K_inj=0,1_Tint=1250_low=0,70_tdead=0,25_dt=1e-4/coherence_resonance_beta={:f}".format(b) + ".txt", skiprows=2)
		bins = np.logspace(np.log10(binsLowLim), np.log10(binsUpLim), 100)
		plt.hist(TISI, bins=bins, density=True, log=True, alpha=0.85, label=r"$\beta=$" + str(b), color=colorList[i])
		i += 1
	
	plt.xscale("log")
	plt.xlabel(r"interspike interval time $T_{ISI}$ / ns")
	plt.ylabel(r"frequency of $T_{ISI,i}$", labelpad=-1.0)
	plt.yticks([0.0001, 0.01, 1.0])
	plt.xlim(binsLowLim, binsUpLim)
	plt.legend(loc="upper right")
	plt.show()


# ~ plotCoherenceResonanceInjection()
# ~ plotTISIHistogramInjection(0.100000, 0.010000, 0.001000)
# ~ plotCoherenceResonanceInjectionFeedback()
# ~ plotTISIHistogramInjectionFeedback(0.523299, 0.052330, 0.005233)
plotCoherenceResonanceNetwork()

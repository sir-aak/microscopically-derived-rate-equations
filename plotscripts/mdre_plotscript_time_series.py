import numpy as np
import matplotlib.pyplot as plt


# plots time series for MDRE model
def plotTimeSeries ():
	
	with open("time_series.txt", "r") as file:
		lines = file.readlines()
	
	time                = []
	intensity           = []
	E_real              = []
	E_imag              = []
	amplitude           = []
	phase               = []
	rho_GS_e_act        = []
	rho_GS_h_act        = []
	rho_GS_e_inact      = []
	rho_GS_h_inact      = []
	rho_ES_e            = []
	rho_ES_h            = []
	w_e                 = []
	w_h                 = []
	electronSum         = []
	holeSum             = []
	chargeConservation  = []
	dchargeConservation = []
	occupation_GS       = []
	occupation_ES       = []
	
	for line in lines:
		time.append(float((line.split('	')[0])))
		intensity.append(float((line.split('	')[1])))
		E_real.append(float((line.split('	')[2])))
		E_imag.append(float((line.split('	')[3])))
		amplitude.append(float((line.split('	')[4])))
		phase.append(float((line.split('	')[5])))
		rho_GS_e_act.append(float((line.split('	')[6])))
		rho_GS_h_act.append(float((line.split('	')[7])))
		rho_GS_e_inact.append(float((line.split('	')[8])))
		rho_GS_h_inact.append(float((line.split('	')[9])))
		rho_ES_e.append(float((line.split('	')[10])))
		rho_ES_h.append(float((line.split('	')[11])))
		w_e.append(float((line.split('	')[12])))
		w_h.append(float((line.split('	')[13])))
		electronSum.append(float((line.split('	')[14])))
		holeSum.append(float((line.split('	')[15])))
		chargeConservation.append(float((line.split('	')[16])))
		dchargeConservation.append(float((line.split('	')[17])))
		occupation_GS.append(float((line.split('	')[18])))
		occupation_ES.append(float((line.split('	')[19])))
	
	time                = np.array(time)
	intensity           = np.array(intensity)
	E_real              = np.array(E_real)
	E_imag              = np.array(E_imag)
	amplitude           = np.array(amplitude)
	phase               = np.array(phase)
	rho_GS_e_act        = np.array(rho_GS_e_act)
	rho_GS_h_act        = np.array(rho_GS_h_act)
	rho_GS_e_inact      = np.array(rho_GS_e_inact)
	rho_GS_h_inact      = np.array(rho_GS_h_inact)
	rho_ES_e            = np.array(rho_ES_e)
	rho_ES_h            = np.array(rho_ES_h)
	w_e                 = np.array(w_e)
	w_h                 = np.array(w_h)
	electronSum         = np.array(electronSum)
	holeSum             = np.array(holeSum)
	chargeConservation  = np.array(chargeConservation)
	dchargeConservation = np.array(dchargeConservation)
	occupation_GS       = np.array(occupation_GS)
	occupation_ES       = np.array(occupation_ES)
	
	
	plt.figure(figsize=(10, 10))
	plt.suptitle("MDRE time series", fontsize=14)
	plt.subplots_adjust(hspace=0.50, wspace=0.30)
	
	plt.subplot(7, 2, 1)
	plt.plot(time, intensity)
	plt.ylabel(r"$|E|^2$")
	plt.grid(color="lightgray")
	
	plt.subplot(7, 2, 2)
	plt.plot(time, amplitude)
	plt.ylabel(r"$|E|$")
	plt.grid(color="lightgray")
	
	plt.subplot(7, 2, 3)
	plt.plot(time, electronSum)
	plt.ylabel(r"$\Sigma e^-$")
	plt.grid(color="lightgray")
	
	plt.subplot(7, 2, 4)
	plt.plot(time, holeSum)
	plt.ylabel(r"$\Sigma h$")
	plt.grid(color="lightgray")
	
	plt.subplot(7, 2, 5)
	plt.plot(time, chargeConservation)
	plt.ylabel(r"$\frac{\Sigma e^- - \Sigma h}{\Sigma e^-}$")
	plt.grid(color="lightgray")
	
	plt.subplot(7, 2, 6)
	plt.plot(time, dchargeConservation)
	plt.ylabel(r"$\frac{1}{\Sigma e^-}$ $\frac{d}{dt}$ ($\Sigma e^- - \Sigma h$) / ns$^{-1}$")
	plt.grid(color="lightgray")
	
	plt.subplot(7, 2, 7)
	plt.plot(time, rho_GS_e_act)
	plt.ylabel(r"$\rho_{GS,e}^{act}$")
	plt.ylim(0.0, 1.0)
	plt.grid(color="lightgray")
	
	plt.subplot(7, 2, 8)
	plt.plot(time, rho_GS_h_act)
	plt.ylabel(r"$\rho_{GS,h}^{act}$")
	plt.ylim(0.0, 1.0)
	plt.grid(color="lightgray")
	
	plt.subplot(7, 2, 9)
	plt.plot(time, rho_GS_e_inact)
	plt.ylabel(r"$\rho_{GS,e}^{inact}$")
	plt.ylim(0.0, 1.0)
	plt.grid(color="lightgray")
	
	plt.subplot(7, 2, 10)
	plt.plot(time, rho_GS_h_inact)
	plt.ylabel(r"$\rho_{GS,h}^{inact}$")
	plt.ylim(0.0, 1.0)
	plt.grid(color="lightgray")
	
	plt.subplot(7, 2, 11)
	plt.plot(time, rho_ES_e)
	plt.ylabel(r"$\rho_{ES,e}$")
	plt.ylim(0.0, 1.0)
	plt.grid(color="lightgray")
	
	plt.subplot(7, 2, 12)
	plt.plot(time, rho_ES_h)
	plt.ylabel(r"$\rho_{ES,h}$")
	plt.ylim(0.0, 1.0)
	plt.grid(color="lightgray")
	
	plt.subplot(7, 2, 13)
	plt.plot(time, w_e)
	plt.xlabel(r"time $t$ / ns")
	plt.ylabel(r"$w_e$")
	plt.grid(color="lightgray")
	
	plt.subplot(7, 2, 14)
	plt.plot(time, w_h)
	plt.xlabel(r"time $t$ / ns")
	plt.ylabel(r"$w_h$")
	plt.grid(color="lightgray")
	
	plt.show()


plotTimeSeries()


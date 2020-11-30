import numpy as np
import matplotlib.pyplot as plt


# plots time series for MDRE model 
# comparison between two numerical integration methods
def plotTimeSeries ():
	
	with open("time_series_mdre_method1.txt", "r") as file_method1:
		Lines_method1 = file_method1.readlines()
	
	time_method1                = []
	intensity_method1           = []
	amplitude_method1           = []
	phase_method1               = []
	rho_GS_e_act_method1        = []
	rho_GS_h_act_method1        = []
	rho_GS_e_inact_method1      = []
	rho_GS_h_inact_method1      = []
	rho_ES_e_method1            = []
	rho_ES_h_method1            = []
	w_e_method1                 = []
	w_h_method1                 = []
	electronSum_method1         = []
	holeSum_method1             = []
	chargeConservation_method1  = []
	dchargeConservation_method1 = []
	
	for Line_method1 in Lines_method1:
		
		time_method1.append(float((Line_method1.split('	')[0])))
		intensity_method1.append(float((Line_method1.split('	')[1])))
		amplitude_method1.append(float((Line_method1.split('	')[2])))
		phase_method1.append(float((Line_method1.split('	')[3])))
		rho_GS_e_act_method1.append(float((Line_method1.split('	')[4])))
		rho_GS_h_act_method1.append(float((Line_method1.split('	')[5])))
		rho_GS_e_inact_method1.append(float((Line_method1.split('	')[6])))
		rho_GS_h_inact_method1.append(float((Line_method1.split('	')[7])))
		rho_ES_e_method1.append(float((Line_method1.split('	')[8])))
		rho_ES_h_method1.append(float((Line_method1.split('	')[9])))
		w_e_method1.append(float((Line_method1.split('	')[10])))
		w_h_method1.append(float((Line_method1.split('	')[11])))
		electronSum_method1.append(float((Line_method1.split('	')[12])))
		holeSum_method1.append(float((Line_method1.split('	')[13])))
		chargeConservation_method1.append(float((Line_method1.split('	')[14])))
		dchargeConservation_method1.append(float((Line_method1.split('	')[15])))
	
	time_method1                = np.array(time_method1)
	intensity_method1           = np.array(intensity_method1)
	amplitude_method1           = np.array(amplitude_method1)
	phase_method1               = np.array(phase_method1)
	rho_GS_e_act_method1        = np.array(rho_GS_e_act_method1)
	rho_GS_h_act_method1        = np.array(rho_GS_h_act_method1)
	rho_GS_e_inact_method1      = np.array(rho_GS_e_inact_method1)
	rho_GS_h_inact_method1      = np.array(rho_GS_h_inact_method1)
	rho_ES_e_method1            = np.array(rho_ES_e_method1)
	rho_ES_h_method1            = np.array(rho_ES_h_method1)
	w_e_method1                 = np.array(w_e_method1)
	w_h_method1                 = np.array(w_h_method1)
	electronSum_method1         = np.array(electronSum_method1)
	holeSum_method1             = np.array(holeSum_method1)
	chargeConservation_method1  = np.array(chargeConservation_method1)
	dchargeConservation_method1 = np.array(dchargeConservation_method1)
	
	
	with open("time_series_mdre_method2.txt", "r") as file_method2:
		Lines_method2 = file_method2.readlines()
	
	time_method2                = []
	intensity_method2           = []
	amplitude_method2           = []
	phase_method2               = []
	rho_GS_e_act_method2        = []
	rho_GS_h_act_method2        = []
	rho_GS_e_inact_method2      = []
	rho_GS_h_inact_method2      = []
	rho_ES_e_method2            = []
	rho_ES_h_method2            = []
	w_e_method2                 = []
	w_h_method2                 = []
	electronSum_method2         = []
	holeSum_method2             = []
	chargeConservation_method2  = []
	dchargeConservation_method2 = []
	
	for Line_method2 in Lines_method2:
		
		time_method2.append(float((Line_method2.split('	')[0])))
		intensity_method2.append(float((Line_method2.split('	')[1])))
		amplitude_method2.append(float((Line_method2.split('	')[2])))
		phase_method2.append(float((Line_method2.split('	')[3])))
		rho_GS_e_act_method2.append(float((Line_method2.split('	')[4])))
		rho_GS_h_act_method2.append(float((Line_method2.split('	')[5])))
		rho_GS_e_inact_method2.append(float((Line_method2.split('	')[6])))
		rho_GS_h_inact_method2.append(float((Line_method2.split('	')[7])))
		rho_ES_e_method2.append(float((Line_method2.split('	')[8])))
		rho_ES_h_method2.append(float((Line_method2.split('	')[9])))
		w_e_method2.append(float((Line_method2.split('	')[10])))
		w_h_method2.append(float((Line_method2.split('	')[11])))
		electronSum_method2.append(float((Line_method2.split('	')[12])))
		holeSum_method2.append(float((Line_method2.split('	')[13])))
		chargeConservation_method2.append(float((Line_method2.split('	')[14])))
		dchargeConservation_method2.append(float((Line_method2.split('	')[15])))
	
	time_method2                = np.array(time_method2)
	intensity_method2           = np.array(intensity_method2)
	amplitude_method2           = np.array(amplitude_method2)
	phase_method2               = np.array(phase_method2)
	rho_GS_e_act_method2        = np.array(rho_GS_e_act_method2)
	rho_GS_h_act_method2        = np.array(rho_GS_h_act_method2)
	rho_GS_e_inact_method2      = np.array(rho_GS_e_inact_method2)
	rho_GS_h_inact_method2      = np.array(rho_GS_h_inact_method2)
	rho_ES_e_method2            = np.array(rho_ES_e_method2)
	rho_ES_h_method2            = np.array(rho_ES_h_method2)
	w_e_method2                 = np.array(w_e_method2)
	w_h_method2                 = np.array(w_h_method2)
	electronSum_method2         = np.array(electronSum_method2)
	holeSum_method2             = np.array(holeSum_method2)
	chargeConservation_method2  = np.array(chargeConservation_method2)
	dchargeConservation_method2 = np.array(dchargeConservation_method2)
	
	
	plt.figure
	plt.suptitle(r"MDRE time series: implicit Euler ($dt=10^{-5}$ ns) vs. Runge-Kutta 4 ($dt=10^{-3}$ ns)")
	plt.subplots_adjust(hspace=0.50, wspace=0.30)
	
	plt.subplot(7, 2, 1)
	plt.plot(time_method1, intensity_method1, color="blue", alpha=0.5, label="implicit Euler")
	plt.plot(time_method2, intensity_method2, color="red", alpha=0.5, label="Runge-Kutta 4")
	plt.ylabel(r"$|E|^2$ / cm$^{-2}$")
	plt.grid(color="lightgray")
	plt.legend(bbox_to_anchor=(1, 1), loc="upper right")
	
	plt.subplot(7, 2, 2)
	plt.plot(time_method1, amplitude_method1, color="blue", alpha=0.5)
	plt.plot(time_method2, amplitude_method2, color="red", alpha=0.5)
	plt.ylabel(r"$|E|$ / cm$^{-1}$")
	plt.grid(color="lightgray")
	
	plt.subplot(7, 2, 3)
	plt.plot(time_method1, electronSum_method1, color="blue", alpha=0.5)
	plt.plot(time_method2, electronSum_method2, color="red", alpha=0.5)
	plt.ylabel(r"$\Sigma e^-$ / cm$^{-2}$")
	plt.grid(color="lightgray")
	
	plt.subplot(7, 2, 4)
	plt.plot(time_method1, holeSum_method1, color="blue", alpha=0.5)
	plt.plot(time_method2, holeSum_method2, color="red", alpha=0.5)
	plt.ylabel(r"$\Sigma h$ / cm$^{-2}$")
	plt.grid(color="lightgray")
	
	plt.subplot(7, 2, 5)
	plt.plot(time_method1, chargeConservation_method1, color="blue", alpha=0.5)
	plt.plot(time_method2, chargeConservation_method2, color="red", alpha=0.5)
	plt.ylabel(r"$\frac{\Sigma e^- - \Sigma h}{\Sigma e^-}$")
	plt.grid(color="lightgray")
	
	plt.subplot(7, 2, 6)
	plt.plot(time_method1, dchargeConservation_method1, color="blue", alpha=0.5)
	plt.plot(time_method2, dchargeConservation_method2, color="red", alpha=0.5)
	plt.ylabel(r"$\frac{1}{\Sigma e^-}$ $\frac{d}{dt}$ ($\Sigma e^- - \Sigma h$) / ns$^{-1}$")
	plt.grid(color="lightgray")
	
	plt.subplot(7, 2, 7)
	plt.plot(time_method1, rho_GS_e_act_method1, color="blue", alpha=0.5)
	plt.plot(time_method2, rho_GS_e_act_method2, color="red", alpha=0.5)
	plt.ylabel(r"$\rho_{GS,e}^{act}$")
	plt.ylim(0.0, 1.0)
	plt.grid(color="lightgray")
	
	plt.subplot(7, 2, 8)
	plt.plot(time_method1, rho_GS_h_act_method1, color="blue", alpha=0.5)
	plt.plot(time_method2, rho_GS_h_act_method2, color="red", alpha=0.5)
	plt.ylabel(r"$\rho_{GS,h}^{act}$")
	plt.ylim(0.0, 1.0)
	plt.grid(color="lightgray")
	
	plt.subplot(7, 2, 9)
	plt.plot(time_method1, rho_GS_e_inact_method1, color="blue", alpha=0.5)
	plt.plot(time_method2, rho_GS_e_inact_method2, color="red", alpha=0.5)
	plt.ylabel(r"$\rho_{GS,e}^{inact}$")
	plt.ylim(0.0, 1.0)
	plt.grid(color="lightgray")
	
	plt.subplot(7, 2, 10)
	plt.plot(time_method1, rho_GS_h_inact_method1, color="blue", alpha=0.5)
	plt.plot(time_method2, rho_GS_h_inact_method2, color="red", alpha=0.5)
	plt.ylabel(r"$\rho_{GS,h}^{inact}$")
	plt.ylim(0.0, 1.0)
	plt.grid(color="lightgray")
	
	plt.subplot(7, 2, 11)
	plt.plot(time_method1, rho_ES_e_method1, color="blue", alpha=0.5)
	plt.plot(time_method2, rho_ES_e_method2, color="red", alpha=0.5)
	plt.ylabel(r"$\rho_{ES,e}$")
	plt.ylim(0.0, 1.0)
	plt.grid(color="lightgray")
	
	plt.subplot(7, 2, 12)
	plt.plot(time_method1, rho_ES_h_method1, color="blue", alpha=0.5)
	plt.plot(time_method2, rho_ES_h_method2, color="red", alpha=0.5)
	plt.ylabel(r"$\rho_{ES,h}$")
	plt.ylim(0.0, 1.0)
	plt.grid(color="lightgray")
	
	plt.subplot(7, 2, 13)
	plt.plot(time_method1, w_e_method1, color="blue", alpha=0.5)
	plt.plot(time_method2, w_e_method2, color="red", alpha=0.5)
	plt.xlabel(r"time $t$ / ns")
	plt.ylabel(r"$w_e$ / cm$^{-2}$")
	plt.grid(color="lightgray")
	
	plt.subplot(7, 2, 14)
	plt.plot(time_method1, w_h_method1, color="blue", alpha=0.5)
	plt.plot(time_method2, w_h_method2, color="red", alpha=0.5)
	plt.xlabel(r"time $t$ / ns")
	plt.ylabel(r"$w_h$ / cm$^{-2}$")
	plt.grid(color="lightgray")
	
	plt.show()


plotTimeSeries()

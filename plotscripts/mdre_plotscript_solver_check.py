import numpy as np
import matplotlib.pyplot as plt


# plots time series for MDRE model 
# comparison between Andrej's and Stefan's time series
def plotTimeSeries ():
	
	# extraction of Andrej's data from file
	with open("time_series_mdre_andrej.txt", "r") as file_andrej:
		Lines_andrej = file_andrej.readlines()
	
	time_andrej                = []
	intensity_andrej           = []
	amplitude_andrej           = []
	phase_andrej               = []
	rho_GS_e_act_andrej        = []
	rho_GS_h_act_andrej        = []
	rho_GS_e_inact_andrej      = []
	rho_GS_h_inact_andrej      = []
	rho_ES_e_andrej            = []
	rho_ES_h_andrej            = []
	w_e_andrej                 = []
	w_h_andrej                 = []
	
	for Line_andrej in Lines_andrej:
		
		time_andrej.append(float((Line_andrej.split('	')[0])))
		intensity_andrej.append(float((Line_andrej.split('	')[1])))
		amplitude_andrej.append(float((Line_andrej.split('	')[2])))
		phase_andrej.append(float((Line_andrej.split('	')[3])))
		rho_GS_e_act_andrej.append(float((Line_andrej.split('	')[4])))
		rho_GS_h_act_andrej.append(float((Line_andrej.split('	')[5])))
		rho_GS_e_inact_andrej.append(float((Line_andrej.split('	')[6])))
		rho_GS_h_inact_andrej.append(float((Line_andrej.split('	')[7])))
		rho_ES_e_andrej.append(float((Line_andrej.split('	')[8])))
		rho_ES_h_andrej.append(float((Line_andrej.split('	')[9])))
		w_e_andrej.append(float((Line_andrej.split('	')[10])))
		w_h_andrej.append(float((Line_andrej.split('	')[11])))
	
	time_andrej                = np.array(time_andrej)
	intensity_andrej           = np.array(intensity_andrej)
	amplitude_andrej           = np.array(amplitude_andrej)
	phase_andrej               = np.array(phase_andrej)
	rho_GS_e_act_andrej        = np.array(rho_GS_e_act_andrej)
	rho_GS_h_act_andrej        = np.array(rho_GS_h_act_andrej)
	rho_GS_e_inact_andrej      = np.array(rho_GS_e_inact_andrej)
	rho_GS_h_inact_andrej      = np.array(rho_GS_h_inact_andrej)
	rho_ES_e_andrej            = np.array(rho_ES_e_andrej)
	rho_ES_h_andrej            = np.array(rho_ES_h_andrej)
	w_e_andrej                 = np.array(w_e_andrej)
	w_h_andrej                 = np.array(w_h_andrej)
	
	
	# extraction of Stefan's data from file
	with open("time_series_mdre_stefan.txt", "r") as file_stefan:
		Lines_stefan = file_stefan.readlines()
	
	time_stefan                = []
	amplitude_stefan           = []
	rho_GS_e_act_stefan        = []
	rho_GS_h_act_stefan        = []
	rho_GS_e_inact_stefan      = []
	rho_GS_h_inact_stefan      = []
	rho_ES_e_stefan            = []
	rho_ES_h_stefan            = []
	w_e_stefan                 = []
	w_h_stefan                 = []
	
	for Line_stefan in Lines_stefan[1:]:
		
		time_stefan.append(float((Line_stefan.split('	')[0])))
		amplitude_stefan.append(float((Line_stefan.split('	')[1])))
		rho_GS_e_act_stefan.append(float((Line_stefan.split('	')[3])))
		rho_GS_h_act_stefan.append(float((Line_stefan.split('	')[4])))
		rho_GS_e_inact_stefan.append(float((Line_stefan.split('	')[7])))
		rho_GS_h_inact_stefan.append(float((Line_stefan.split('	')[8])))
		rho_ES_e_stefan.append(float((Line_stefan.split('	')[5])))
		rho_ES_h_stefan.append(float((Line_stefan.split('	')[6])))
		w_e_stefan.append(float((Line_stefan.split('	')[11])))
		w_h_stefan.append(float((Line_stefan.split('	')[12])))
	
	time_stefan                = np.array(time_stefan)
	amplitude_stefan           = np.array(amplitude_stefan)
	rho_GS_e_act_stefan        = np.array(rho_GS_e_act_stefan)
	rho_GS_h_act_stefan        = np.array(rho_GS_h_act_stefan)
	rho_GS_e_inact_stefan      = np.array(rho_GS_e_inact_stefan)
	rho_GS_h_inact_stefan      = np.array(rho_GS_h_inact_stefan)
	rho_ES_e_stefan            = np.array(rho_ES_e_stefan)
	rho_ES_h_stefan            = np.array(rho_ES_h_stefan)
	w_e_stefan                 = np.array(w_e_stefan)
	w_h_stefan                 = np.array(w_h_stefan)
	
	intensity_stefan = np.square(amplitude_stefan)
	nSteps           = intensity_stefan.size
	
	print("\n")
	
	# differences in time series calculated with Euclidean norm
	print("normalized differences Stefan - Andrej, Euclidean norm:")
	print("intensity   : " + str(np.linalg.norm(intensity_stefan - intensity_andrej) / nSteps))
	print("amplitude   : " + str(np.linalg.norm(amplitude_stefan - amplitude_andrej) / nSteps))
	print("rhoGSeact   : " + str(np.linalg.norm(rho_GS_e_act_stefan - rho_GS_e_act_andrej) / nSteps))
	print("rhoGShact   : " + str(np.linalg.norm(rho_GS_h_act_stefan - rho_GS_h_act_andrej) / nSteps))
	print("rhoGSeinact : " + str(np.linalg.norm(rho_GS_e_inact_stefan - rho_GS_e_inact_andrej) / nSteps))
	print("rhoGShinact : " + str(np.linalg.norm(rho_GS_h_inact_stefan - rho_GS_h_inact_andrej) / nSteps))
	print("rhoESe      : " + str(np.linalg.norm(rho_ES_e_stefan - rho_ES_e_andrej) / nSteps))
	print("rhoESh      : " + str(np.linalg.norm(rho_ES_h_stefan - rho_ES_h_andrej) / nSteps))
	print("we          : " + str(np.linalg.norm(w_e_stefan - w_e_andrej) / nSteps))
	print("wh          : " + str(np.linalg.norm(w_h_stefan - w_h_andrej) / nSteps))
	
	print("\n")
	
	# differences in time series calculated with maximum-norm
	print("differences Stefan - Andrej, maximum-norm:")
	print("intensity   : " + str(np.linalg.norm(intensity_stefan - intensity_andrej, np.inf)))
	print("amplitude   : " + str(np.linalg.norm(amplitude_stefan - amplitude_andrej, np.inf)))
	print("rhoGSeact   : " + str(np.linalg.norm(rho_GS_e_act_stefan - rho_GS_e_act_andrej, np.inf)))
	print("rhoGShact   : " + str(np.linalg.norm(rho_GS_h_act_stefan - rho_GS_h_act_andrej, np.inf)))
	print("rhoGSeinact : " + str(np.linalg.norm(rho_GS_e_inact_stefan - rho_GS_e_inact_andrej, np.inf)))
	print("rhoGShinact : " + str(np.linalg.norm(rho_GS_h_inact_stefan - rho_GS_h_inact_andrej, np.inf)))
	print("rhoESe      : " + str(np.linalg.norm(rho_ES_e_stefan - rho_ES_e_andrej, np.inf)))
	print("rhoESh      : " + str(np.linalg.norm(rho_ES_h_stefan - rho_ES_h_andrej, np.inf)))
	print("we          : " + str(np.linalg.norm(w_e_stefan - w_e_andrej, np.inf)))
	print("wh          : " + str(np.linalg.norm(w_h_stefan - w_h_andrej, np.inf)))
	
	print("\n")
	
	
	plt.figure()
	plt.suptitle(r"MDRE time series: Andrejs Runge-Kutta 4 ($dt=10^{-3}$ ns) vs. Stefans Runge-Kutta 4 ($dt=10^{-3}$ ns)")
	plt.subplots_adjust(hspace=0.2, wspace=0.30)
	
	plt.subplot(5, 2, 1)
	plt.plot(time_andrej, intensity_andrej, color="blue", alpha=0.5, label="Andrej")
	plt.plot(time_stefan, intensity_stefan, color="red", alpha=0.5, label="Stefan")
	plt.ylabel(r"$|E|^2$")
	plt.grid(color="lightgray")
	plt.legend(bbox_to_anchor=(1, 1), loc="upper right")
	
	plt.subplot(5, 2, 2)
	plt.plot(time_andrej, amplitude_andrej, color="blue", alpha=0.5)
	plt.plot(time_stefan, amplitude_stefan, color="red", alpha=0.5)
	plt.ylabel(r"$|E|$")
	plt.grid(color="lightgray")
	
	plt.subplot(5, 2, 3)
	plt.plot(time_andrej, rho_GS_e_act_andrej, color="blue", alpha=0.5)
	plt.plot(time_stefan, rho_GS_e_act_stefan, color="red", alpha=0.5)
	plt.ylabel(r"$\rho_{GS,e}^{act}$")
	plt.ylim(0.0, 1.0)
	plt.grid(color="lightgray")
	
	plt.subplot(5, 2, 4)
	plt.plot(time_andrej, rho_GS_h_act_andrej, color="blue", alpha=0.5)
	plt.plot(time_stefan, rho_GS_h_act_stefan, color="red", alpha=0.5)
	plt.ylabel(r"$\rho_{GS,h}^{act}$")
	plt.ylim(0.0, 1.0)
	plt.grid(color="lightgray")
	
	plt.subplot(5, 2, 5)
	plt.plot(time_andrej, rho_GS_e_inact_andrej, color="blue", alpha=0.5)
	plt.plot(time_stefan, rho_GS_e_inact_stefan, color="red", alpha=0.5)
	plt.ylabel(r"$\rho_{GS,e}^{inact}$")
	plt.ylim(0.0, 1.0)
	plt.grid(color="lightgray")
	
	plt.subplot(5, 2, 6)
	plt.plot(time_andrej, rho_GS_h_inact_andrej, color="blue", alpha=0.5)
	plt.plot(time_stefan, rho_GS_h_inact_stefan, color="red", alpha=0.5)
	plt.ylabel(r"$\rho_{GS,h}^{inact}$")
	plt.ylim(0.0, 1.0)
	plt.grid(color="lightgray")
	
	plt.subplot(5, 2, 7)
	plt.plot(time_andrej, rho_ES_e_andrej, color="blue", alpha=0.5)
	plt.plot(time_stefan, rho_ES_e_stefan, color="red", alpha=0.5)
	plt.ylabel(r"$\rho_{ES,e}$")
	plt.ylim(0.0, 1.0)
	plt.grid(color="lightgray")
	
	plt.subplot(5, 2, 8)
	plt.plot(time_andrej, rho_ES_h_andrej, color="blue", alpha=0.5)
	plt.plot(time_stefan, rho_ES_h_stefan, color="red", alpha=0.5)
	plt.ylabel(r"$\rho_{ES,h}$")
	plt.ylim(0.0, 1.0)
	plt.grid(color="lightgray")
	
	plt.subplot(5, 2, 9)
	plt.plot(time_andrej, w_e_andrej, color="blue", alpha=0.5)
	plt.plot(time_stefan, w_e_stefan, color="red", alpha=0.5)
	plt.xlabel(r"time $t$ / ns")
	plt.ylabel(r"$w_e$")
	plt.grid(color="lightgray")
	
	plt.subplot(5, 2, 10)
	plt.plot(time_andrej, w_h_andrej, color="blue", alpha=0.5)
	plt.plot(time_stefan, w_h_stefan, color="red", alpha=0.5)
	plt.xlabel(r"time $t$ / ns")
	plt.ylabel(r"$w_h$")
	plt.grid(color="lightgray")
	
	plt.show()
	
	plt.figure()
	plt.semilogy(time_andrej, amplitude_andrej, color="blue", alpha=0.5, label="Andrej")
	plt.semilogy(time_stefan, amplitude_stefan, color="red", alpha=0.5, label="Stefan")
	plt.xlabel(r"time $t$ / ns")
	plt.ylabel(r"$|E|$")
	plt.grid(color="lightgray")
	plt.legend(loc="lower right")
	
	plt.show()


plotTimeSeries()

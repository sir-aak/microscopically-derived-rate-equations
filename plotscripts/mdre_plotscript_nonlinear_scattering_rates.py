import numpy as np
import re
import matplotlib.pyplot as plt


# plots nonlinear scattering rates in dependence on 2D carrier-density for MDRE model
def plotNonlinearScatteringRates ():
	
	with open("nonlinear_scattering_rates_mdre.txt", "r") as file:
		lines = file.readlines()
	
	w              = []
	S_GS_e_cap_in  = []
	S_GS_h_cap_in  = []
	S_ES_e_cap_in  = []
	S_ES_h_cap_in  = []
	S_e_rel_in     = []
	S_h_rel_in     = []
	S_GS_e_cap_out = []
	S_GS_h_cap_out = []
	S_ES_e_cap_out = []
	S_ES_h_cap_out = []
	S_e_rel_out    = []
	S_h_rel_out    = []
	
	for line in lines[:-1]:			# last line is superfluous
		
		parsedLine = re.sub(" +", "	", line.lstrip())
		
		w.append(float(parsedLine.split('	')[0]))
		S_GS_e_cap_in.append(float(parsedLine.split('	')[1]))
		S_GS_h_cap_in.append(float(parsedLine.split('	')[2]))
		S_ES_e_cap_in.append(float(parsedLine.split('	')[3]))
		S_ES_h_cap_in.append(float(parsedLine.split('	')[4]))
		S_e_rel_in.append(float(parsedLine.split('	')[5]))
		S_h_rel_in.append(float(parsedLine.split('	')[6]))
		S_GS_e_cap_out.append(float(parsedLine.split('	')[7]))
		S_GS_h_cap_out.append(float(parsedLine.split('	')[8]))
		S_ES_e_cap_out.append(float(parsedLine.split('	')[9]))
		S_ES_h_cap_out.append(float(parsedLine.split('	')[10]))
		S_e_rel_out.append(float(parsedLine.split('	')[11]))
		S_h_rel_out.append(float(parsedLine.split('	')[12]))
	
	w              = np.array(w)
	S_GS_e_cap_in  = np.array(S_GS_e_cap_in)
	S_GS_h_cap_in  = np.array(S_GS_h_cap_in)
	S_ES_e_cap_in  = np.array(S_ES_e_cap_in)
	S_ES_h_cap_in  = np.array(S_ES_h_cap_in)
	S_e_rel_in     = np.array(S_e_rel_in)
	S_h_rel_in     = np.array(S_h_rel_in)
	S_GS_e_cap_out = np.array(S_GS_e_cap_out)
	S_GS_h_cap_out = np.array(S_GS_h_cap_out)
	S_ES_e_cap_out = np.array(S_ES_e_cap_out)
	S_ES_h_cap_out = np.array(S_ES_h_cap_out)
	S_e_rel_out    = np.array(S_e_rel_out)
	S_h_rel_out    = np.array(S_h_rel_out)
	
	plt.figure
	plt.suptitle("MDRE nonlinear scattering rates / ns$^{-1}$")
	plt.subplots_adjust(hspace=0.50, wspace=0.30)
	
	plt.subplot(6, 2, 1)
	plt.plot(w, S_GS_e_cap_in)
	plt.ylabel(r"$S_{GS,e}^{cap,in}$")
	plt.grid(color="lightgray")
	
	plt.subplot(6, 2, 2)
	plt.plot(w, S_GS_h_cap_in)
	plt.ylabel(r"$S_{GS,h}^{cap,in}$")
	plt.grid(color="lightgray")
	
	plt.subplot(6, 2, 3)
	plt.plot(w, S_ES_e_cap_in)
	plt.ylabel(r"$S_{ES,e}^{cap,in}$")
	plt.grid(color="lightgray")
	
	plt.subplot(6, 2, 4)
	plt.plot(w, S_ES_h_cap_in)
	plt.ylabel(r"$S_{ES,h}^{cap,in}$")
	plt.grid(color="lightgray")
	
	plt.subplot(6, 2, 5)
	plt.plot(w, S_e_rel_in)
	plt.ylabel(r"$S_e^{rel,in}$")
	plt.grid(color="lightgray")
	
	plt.subplot(6, 2, 6)
	plt.plot(w, S_h_rel_in)
	plt.ylabel(r"$S_h^{rel,in}$")
	plt.grid(color="lightgray")
	
	plt.subplot(6, 2, 7)
	plt.plot(w, S_GS_e_cap_out)
	plt.ylabel(r"$S_{GS,e}^{cap,out}$")
	plt.grid(color="lightgray")
	
	plt.subplot(6, 2, 8)
	plt.plot(w, S_GS_h_cap_out)
	plt.ylabel(r"$S_{GS,h}^{cap,out}$")
	plt.grid(color="lightgray")
	
	plt.subplot(6, 2, 9)
	plt.plot(w, S_ES_e_cap_out)
	plt.ylabel(r"$S_{ES,e}^{cap,out}$")
	plt.grid(color="lightgray")
	
	plt.subplot(6, 2, 10)
	plt.plot(w, S_ES_h_cap_out)
	plt.ylabel(r"$S_{ES,h}^{cap,out}$")
	plt.grid(color="lightgray")
	
	plt.subplot(6, 2, 11)
	plt.plot(w, S_e_rel_out)
	plt.xlabel(r"$w_e$")
	plt.ylabel(r"$S_e^{rel,out}$")
	plt.grid(color="lightgray")
	
	plt.subplot(6, 2, 12)
	plt.plot(w, S_h_rel_out)
	plt.xlabel(r"$w_h$")
	plt.ylabel(r"$S_h^{rel,out}$")
	plt.grid(color="lightgray")
	
	plt.show()


plotNonlinearScatteringRates()

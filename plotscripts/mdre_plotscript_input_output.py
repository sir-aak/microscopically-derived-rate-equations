import numpy as np
import re
import matplotlib.pyplot as plt


# plots input-output curves for MDRE model
def plotInputOutput ():
	
	# read data from file
	with open("input_output_mdre_new.txt", "r") as file:
		lines = file.readlines()
	
	# declaration of empty lists
	current        = []
	intensity      = []
	rho_GS_e_act   = []
	rho_GS_h_act   = []
	rho_GS_e_inact = []
	rho_GS_h_inact = []
	rho_ES_e       = []
	rho_ES_h       = []
	w_e            = []
	w_h            = []
	
	# fill empty lists with data columns
	for line in lines[:-1]:
		
		parsedLine = re.sub(" +", "	", line.lstrip())
		
		current.append(float((parsedLine.split('	')[0])))
		intensity.append(float((parsedLine.split('	')[1])))
		rho_GS_e_act.append(float((parsedLine.split('	')[2])))
		rho_GS_h_act.append(float((parsedLine.split('	')[3])))
		rho_GS_e_inact.append(float((parsedLine.split('	')[4])))
		rho_GS_h_inact.append(float((parsedLine.split('	')[5])))
		rho_ES_e.append(float((parsedLine.split('	')[6])))
		rho_ES_h.append(float((parsedLine.split('	')[7])))
		w_e.append(float((parsedLine.split('	')[8])))
		w_h.append(float((parsedLine.split('	')[9])))
	
	# cast array from filled lists
	current        = np.array(current)
	intensity      = np.array(intensity)
	rho_GS_e_act   = np.array(rho_GS_e_act)
	rho_GS_h_act   = np.array(rho_GS_h_act)
	rho_GS_e_inact = np.array(rho_GS_e_inact)
	rho_GS_h_inact = np.array(rho_GS_h_inact)
	rho_ES_e       = np.array(rho_ES_e)
	rho_ES_h       = np.array(rho_ES_h)
	w_e            = np.array(w_e)
	w_h            = np.array(w_h)
	
	# calculate inversions
	inversion_GS_act   = rho_GS_e_act + rho_GS_h_act - 1.0
	inversion_GS_inact = rho_GS_e_inact + rho_GS_h_inact - 1.0
	inversion_ES       = rho_ES_e + rho_ES_h - 1.0
	
	
	# plot data
	fig, ax1 = plt.subplots()
	# ~ fig.suptitle("MDRE input-output", fontsize=14)
	fig.set_size_inches(5.9, 4.2)
	fig.subplots_adjust(top=0.98, bottom=0.13, left=0.08, right=0.77)
	plt.rcParams.update({"font.size": 10})
	
	#~ ax1.set_title("MDRE input-output")
	# ~ plt.grid(color="lightgray")
	ax2 = ax1.twinx()
	ax3 = ax1.twinx()
	
	ax1.set_xlabel(r"pump current rate $J$ / ns$^{-1}$")
	ax1.set_ylabel(r"steady state intensity $|E|^2$")
	ax1.set_ylim(0.0, 5)
	ax1.set_xlim(-1.0, 51.0)
	
	ax2.set_ylabel(r"steady state 2D charge carrier densities $w_b$", size=10)
	ax2.set_ylim(0.0, 8)
	
	ax3.spines["right"].set_position(("axes", 1.16))
	ax3.set_ylabel(r"steady state inversions $\rho^{(in)act}_{m,e} + \rho^{(in)act}_{m,h}$ - 1", size=10)
	ax3.set_yticks([-1.0, -0.5, 0.0, 0.5, 1.0])
	ax3.set_ylim(-1.0, 1.0)
	
	p1 = ax1.plot(current, intensity, color="crimson", label=r"$|E|^2$")[0]
	p2 = ax3.plot(current, inversion_GS_act, color="orange", label=r"$\rho^{act}_{GS,e} + \rho^{act}_{GS,h} - 1$")[0]
	p3 = ax3.plot(current, inversion_GS_inact, color="gray", linestyle="--", label=r"$\rho^{inact}_{GS,e} + \rho^{inact}_{GS,h} - 1$")[0]
	p4 = ax3.plot(current, inversion_ES, color="cornflowerblue", label=r"$\rho_{ES,e} + \rho_{ES,h} - 1$")[0]
	p5 = ax2.plot(current, w_e, color="green", label=r"$w_e$")[0]
	p6 = ax2.plot(current, w_h, color="green", linestyle="--", label=r"$w_h$")[0]
	
	ax1.yaxis.label.set_color(p1.get_color())
	ax2.yaxis.label.set_color(p5.get_color())
	
	plt.axvline(x=8.5, color="black", linewidth=1.0)
	ax3.text(8.5, 0.675, r"$J_{th}$", ha="center", bbox={"facecolor":"white", "edgecolor":"none", "alpha":1.0})
	
	lines = [p1, p2, p3, p4, p5, p6]
	ax3.legend(lines, [l.get_label() for l in lines], loc="lower right", fancybox=True, framealpha=0.9)
	
	ax1.tick_params(axis="y", colors="crimson")
	ax2.tick_params(axis="y", colors="green")
	
	ax1.set_zorder(1)
	ax1.patch.set_visible(False)
	
	plt.show()


plotInputOutput()


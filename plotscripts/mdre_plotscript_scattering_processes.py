import numpy as np
import re
import matplotlib.pyplot as plt


# plots scattering processes for given initial conditions in dependence on 2D carrier-density for MDRE model
def plotScatteringProcesses ():
    
    with open("scattering_processes_mdre.txt", "r") as file:
        lines = file.readlines()
    
    w                = []
    S_GS_e_cap_act   = []
    S_GS_h_cap_act   = []
    S_GS_e_cap_inact = []
    S_GS_h_cap_inact = []
    S_ES_e_cap       = []
    S_ES_h_cap       = []
    S_e_rel_act      = []
    S_h_rel_act      = []
    S_e_rel_inact    = []
    S_h_rel_inact    = []
    
    for line in lines[:-1]:			# last line is superfluous
        
        parsedLine = re.sub(" +", "	", line.lstrip())
        
        w.append(float(parsedLine.split('	')[0]))
        S_GS_e_cap_act.append(float(parsedLine.split('	')[1]))
        S_GS_h_cap_act.append(float(parsedLine.split('	')[2]))
        S_GS_e_cap_inact.append(float(parsedLine.split('	')[3]))
        S_GS_h_cap_inact.append(float(parsedLine.split('	')[4]))
        S_ES_e_cap.append(float(parsedLine.split('	')[5]))
        S_ES_h_cap.append(float(parsedLine.split('	')[6]))
        S_e_rel_act.append(float(parsedLine.split('	')[7]))
        S_h_rel_act.append(float(parsedLine.split('	')[8]))
        S_e_rel_inact.append(float(parsedLine.split('	')[9]))
        S_h_rel_inact.append(float(parsedLine.split('	')[10]))
    
    w                = np.array(w)
    S_GS_e_cap_act   = np.array(S_GS_e_cap_act)
    S_GS_h_cap_act   = np.array(S_GS_h_cap_act)
    S_GS_e_cap_inact = np.array(S_GS_e_cap_inact)
    S_GS_h_cap_inact = np.array(S_GS_h_cap_inact)
    S_ES_e_cap       = np.array(S_ES_e_cap)
    S_ES_h_cap       = np.array(S_ES_h_cap)
    S_e_rel_act      = np.array(S_e_rel_act)
    S_h_rel_act      = np.array(S_h_rel_act)
    S_e_rel_inact    = np.array(S_e_rel_inact)
    S_h_rel_inact    = np.array(S_h_rel_inact)
    
    plt.figure(figsize=(11.8, 2.5))
#    plt.suptitle("MDRE scattering processes at fixed point / ns$^{-1}$ ")
    plt.subplots_adjust(wspace=0.70, left=0.08, right=0.99, bottom=0.25, top=0.98)
    plt.rcParams.update({'font.size': 18})
    
    plt.subplot(1, 4, 1)
    plt.plot(w, S_GS_e_cap_act, color="orange")
    plt.plot(w, S_GS_h_cap_act, color="orange", linestyle="dashed")
    plt.plot(w, S_GS_e_cap_inact, color="gray")
    plt.plot(w, S_GS_h_cap_inact, color="gray", linestyle="dashed")
    plt.xlabel(r"$w_b$")
    plt.ylabel(r"$S_{GS,b}^{cap,act}$")
    plt.xticks([0.0, 5.0, 10.0])
    
    plt.subplot(1, 4, 2)
    plt.plot(w, S_ES_e_cap, color="cornflowerblue")
    plt.plot(w, S_ES_h_cap, color="cornflowerblue", linestyle="dashed")
    plt.xlabel(r"$w_b$")
    plt.ylabel(r"$S_{ES,b}^{cap}$")
    plt.xticks([0.0, 5.0, 10.0])
    
    plt.subplot(1, 4, 3)
    plt.plot(w, S_e_rel_act, color="orange")
    plt.plot(w, S_h_rel_act, color="orange", linestyle="dashed")
    plt.xlabel(r"$w_b$")
    plt.ylabel(r"$S_b^{rel,act}$")
    plt.xticks([0.0, 5.0, 10.0])
    
    plt.subplot(1, 4, 4)
    plt.plot(w, S_e_rel_inact, color="gray")
    plt.plot(w, S_h_rel_inact, color="gray", linestyle="dashed")
    plt.xlabel(r"$w_b$")
    plt.ylabel(r"$S_b^{rel,inact}$", labelpad=-1.0)
    plt.xticks([0.0, 5.0, 10.0])
    
    plt.show()


plotScatteringProcesses()


import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


# plots intensity fit for MDRE model
def plotIntensityFit ():
	
	with open("time_series_mdre_relaxation_oscillations.txt", "r") as file:
		lines = file.readlines()
	
	time      = []
	intensity = []
	
	for line in lines:
		
		time.append(float((line.split('	')[0])))
		intensity.append(float((line.split('	')[1])))
	
	time      = np.array(time)
	intensity = np.array(intensity)
	
	
	# fit function
	#~ def func (t, a, b, phi, omega, Gamma):
		#~ return (a * np.cos(omega * t + phi) * np.exp(-Gamma * t) + b)
	
	
	def func (t, a, I0, Gamma, omega):
		return (I0 + a * np.exp(-(Gamma + 1j * omega) * t))
	
	
	# initial guess parameters:
	
	# for rates x1 and J = 17 old function
	#~ a     = 4.4e13
	#~ b     = 1.06022
	#~ phi   = 2.0
	#~ omega = 11.66
	#~ Gamma = 10.0
	
	# for rates x0.02 and J = 21 old function
	#~ a     = 2.79e13
	#~ b     = 0.101643
	#~ phi   = 0.86
	#~ omega = 22.76
	#~ Gamma = 20.7
	
	
	# for rates x1 and J = 17 new function
	a     = 6.18e9 + 1j * 8.36e8
	I0    = 1.06
	Gamma = 7.5
	omega = 12.5
	
	
	initial_guess = [a, I0, Gamma, omega]
	
	popt, pcov = curve_fit(func, time, intensity, p0=initial_guess)
	
	#~ np.set_printoptions(precision=2)
	#~ fitted_parameters = "fitted parameters:\n" + "\n$a$ = " + str(popt[0]) + "\n$b$ = " + str(popt[1]) + "\n$\phi$ = " + str(popt[2]) + "\n$\omega$ = " + str(popt[3]) + "\n$\Gamma$ = " + str(popt[4])
	
	#~ plt.figure(figsize=(8, 6))
	#~ plt.suptitle("MDRE function fit for relaxation oscillation parameters", fontsize=14)
	
	#~ plt.plot(time, intensity, color="blue", label="time series")
	#~ plt.plot(time, func(time, *initial_guess), color="gray", linestyle="--", label="guessed function")
	#~ plt.plot(time, func(time, *popt), color="red", label="fitted function")
	
	# for rates x1 and J = 17
	#~ plt.title("rates x1,  $J=17$")
	#~ plt.text(4.125, 2.54, r"fit function: $a$ cos($\omega t + \phi$) exp($-\Gamma t$) + b")
	#~ plt.text(4.125, 1.85, fitted_parameters, wrap=True)
	#~ plt.ylim(0.75, 3.0)
	
	# for rates x 0.02 and J = 21
	#~ plt.title("rates x0.02,  $J=21$")
	#~ plt.text(2.65, 0.195, r"fit function: $a$ cos($\omega t + \phi$) exp($-\Gamma t$) + b")
	#~ plt.text(2.65, 0.155, fitted_parameters, wrap=True)
	#~ plt.ylim(0.08, 0.22)
	
	#~ plt.xlabel("time / ns")
	#~ plt.ylabel("intensity")
	#~ plt.grid(color="lightgray")
	#~ plt.legend(loc="upper right")
	
	#~ plt.show()


plotIntensityFit()

























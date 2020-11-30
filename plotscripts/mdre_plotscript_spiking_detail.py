import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d


# index boundaries for time 3D plot
nStart = 140000
nEnd   = 160000


with open("time_series_stochastic_old.txt", "r") as file:
	lines = file.readlines()

time      = []
intensity = []
E_real    = []
E_imag    = []

for line in lines:
	time.append(float((line.split('	')[0])))
	intensity.append(float((line.split('	')[1])))
	E_real.append(float((line.split('	')[2])))
	E_imag.append(float((line.split('	')[3])))

time      = np.array(time)
intensity = np.array(intensity)
E_real    = np.array(E_real)
E_imag    = np.array(E_imag)


fig, ax = plt.subplots()
fig.set_size_inches(5.9, 4.8)
fig.subplots_adjust(top=0.99, bottom=0.15, left=0.16, right=0.95)
ax.plot(time, intensity, color="darkblue")
ax.set_xlabel(r"time $t$ / ns", fontsize=18.0)
ax.set_ylabel(r"intensity $|E|^2$", fontsize=18.0)
ax.set_xlim(140.0, 160.0)
ax.set_ylim(0.45, 2.05)
ax.set_yticks([0.5, 1.0, 1.5, 2.0])
ax.tick_params(axis="x", labelsize=18.0)
ax.tick_params(axis="y", labelsize=18.0)
ax.grid(color="lightgray")

fig, ax = plt.subplots()
fig.set_size_inches(5.9, 4.8)
plt.rcParams.update({"font.size": 18})
plt.subplots_adjust(top=1.06, bottom=0.05, left=-0.09, right=0.96)
ax = plt.axes(projection="3d")
ax.plot3D(time[nStart:nEnd], E_real[nStart:nEnd], E_imag[nStart:nEnd], color="darkblue")
ax.set_xlabel(r"time $t$ / ns")
ax.set_ylabel(r"Re($E$)")
ax.set_zlabel(r"Im($E$)")
ax.xaxis.labelpad=16
ax.yaxis.labelpad=11
ax.zaxis.labelpad=8
ax.set_xlim(140, 160.0)
ax.set_ylim(-1.5, 1.5)
ax.set_zlim(-1.1, 1.1)
ax.set_xticks([140.0, 145.0, 150.0, 155.0, 160.0])
ax.set_yticks([-1.0, 0.0, 1.0])
ax.set_zticks([-1.0, 0.0, 1.0])

plt.show()


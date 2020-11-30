import numpy as np
# ~ import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d


# ~ mpl.use("TKAGG")

# index boundaries for time 3D plot
nStart = 2500
nEnd   = 5000

data   = np.loadtxt("time_series_rotating_frame.txt") 
time   = data[:, 0]
E_real = data[:, 2]
E_imag = data[:, 3]


fig, ax = plt.subplots()
fig.set_size_inches(5.9, 4.8)
plt.rcParams.update({"font.size": 18})
plt.subplots_adjust(top=1.06, bottom=0.05, left=-0.12, right=0.96)
ax = plt.axes(projection='3d')
ax.plot3D(time[nStart:nEnd], E_real[nStart:nEnd], E_imag[nStart:nEnd], color="black")
ax.set_xlabel(r"time $t$ / ns")
ax.set_ylabel(r"Re($E$)")
ax.set_zlabel(r"Im($E$)")

ax.set_ylim(-1.5, 1.5)
ax.set_zlim(-1.5, 1.5)

ax.xaxis.labelpad=17
ax.yaxis.labelpad=11
ax.zaxis.labelpad=5

ax.set_yticks([-1.0, 0.0, 1.0])
ax.set_zticks([-1.0, 0.0, 1.0])

plt.show()


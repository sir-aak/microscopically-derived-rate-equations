import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors

xbins  = 200
ybins  = 200
bounds = 0.02

data    = np.loadtxt("noise_check.txt", delimiter="	")
dW_real = data[:, 0]
dW_imag = data[:, 1]

numberOfSamples = dW_real.size

plt.figure(figsize=(5.9, 4.5))
plt.subplots_adjust(top=0.98, bottom=0.15, left=0.21, right=0.92)
# ~ plt.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
# ~ plt.suptitle(r"noise check histogram for " + str(numberOfSamples) + " samples", fontsize=14)
# ~ plt.subplots_adjust(wspace=0.30)
# ~ plt.title(r"d$W_t \sim \mathcal{CN}(\mu = 0, \sigma^2 =$d" + r"$t)$", fontsize=18)
plt.hist2d(dW_real, dW_imag, bins=(xbins, ybins), cmap="binary")
plt.xlim(-bounds, bounds)
plt.ylim(-bounds, bounds)
plt.xlabel(r"Re(d$W_t$)", fontsize=18)
plt.ylabel(r"Im(d$W_t$)", fontsize=18, labelpad=-5.0)
plt.xticks([-0.015, 0.0, 0.015], ["-0.015", "0", "0.015"], fontsize=18)
plt.yticks([-0.015, 0.0, 0.015], ["-0.015", "0", "0.015"], fontsize=18)

cb = plt.colorbar()
cb.set_label("frequency", fontsize=18)
cb.ax.tick_params(labelsize=18)

plt.show()



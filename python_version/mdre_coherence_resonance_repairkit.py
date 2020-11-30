import numpy as np
import matplotlib.pyplot as plt

# copy spikeTimes from job file with broken data
# computes mean and variance to repair coherence resonance data
spikeTimes = np.loadtxt("injection/55.txt")
TISI       = np.diff(spikeTimes)

print("mean     = ", np.mean(TISI))
print("variance = ", np.var(TISI))

# ~ plt.plot(TISI)
# ~ plt.show()

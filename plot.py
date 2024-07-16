import numpy as np
import matplotlib.pyplot as plt


plt.rcParams['text.usetex'] = True

bif = np.loadtxt("bif.txt")

plt.figure(figsize = (5, 3), constrained_layout = True)
plt.plot(bif[:, 1], bif[:, 0], ",b")
plt.xlim(0.0, 3.0)
plt.ylim(-2.0, 5.0)
plt.xlabel("$c$", fontsize = 20)
plt.ylabel("$x$", fontsize = 20)
plt.savefig("bif_lorenz.jpg", dpi = 300)

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition

dimfile = "data/dim_2d_lattice_100.dat"
outfile = "plots/out/dim_2d_lattice_100.pdf"

target_dimension = 2
inset_sigma = 1000

print("Reading file {} and writing plot to {}".format(dimfile, outfile))
# 2D grid plot
data = np.loadtxt(dimfile)
startnode, sigma, dim = data.T
startnodes = np.unique(startnode)

fig, ax1 = plt.subplots()

for node in startnodes:
    mask = (node == startnode) & (~np.isnan(dim))
    ax1.plot(sigma[mask], dim[mask], c="tab:blue")

plt.axhline(target_dimension, c="tab:green", ls="--")

plt.xlim([np.min(sigma), np.max(sigma)])
plt.ylim(0, 3)

plt.xlabel("$\\sigma$")
plt.ylabel("$d_{\\rm spec}$")

ax2 = plt.axes([0, 0, 1, 1])
# Manually set the position and relative size of the inset axes within ax1
ip = InsetPosition(ax1, [0.65, 0.65, 0.3, 0.3])
ax2.set_axes_locator(ip)
ax2.set_xlabel("$\\sigma$")
ax2.set_ylabel("$d_{\\rm spec}$")

ax2.plot(sigma[:inset_sigma], dim[:inset_sigma], c="tab:blue")
ax2.axhline(target_dimension, c="tab:green", ls="--")

plt.savefig(outfile)

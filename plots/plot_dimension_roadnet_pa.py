import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition

dimfile = 'data/dim_roadnet_pa.dat'
outfile = 'plots/out/dim_roadnet_pa.pdf'

print("Reading file {} and writing plot to {}".format(dimfile, outfile))
# 1D chain plot
data = np.loadtxt(dimfile)
startnode, sigma, dim = data.T
startnodes = np.unique(startnode)

fig, ax1 = plt.subplots()

for node in startnodes:
    mask = (node == startnode) & (~np.isnan(dim))
    ax1.plot(sigma[mask], dim[mask], alpha=0.2, c='tab:blue')

sigmas = np.unique(sigma)
means = np.zeros(sigmas.shape)
stds = np.zeros(sigmas.shape)

# poor mans implementation of a pivot table
for ind, s in enumerate(sigmas):
    relevant_data = data[data[:, 1] == s]
    dimensions = relevant_data[:, 2]
    means[ind] = np.mean(dimensions)
    stds[ind] = np.std(dimensions)

ax1.plot(sigmas, means, c='tab:orange')
plt.fill_between(sigmas, means + stds, means - stds, alpha=0.2, color='tab:orange')

plt.axhline(2, c='tab:green', ls='--')

plt.xlim([np.min(sigma), np.max(sigma)])
plt.ylim(1, 3.5)

plt.xlabel('$\\sigma$')
plt.ylabel('$d_{\\rm spec}$')

ax2 = plt.axes([0, 0, 1, 1])
# Manually set the position and relative size of the inset axes within ax1
ip = InsetPosition(ax1, [0.65, 0.65, 0.3, 0.3])
ax2.set_axes_locator(ip)
ax2.set_xlabel('$d_{\\rm spec}(\\sigma = 450)$')
ax2.set_ylabel('$n$')
ax2.axvline(np.mean(dimensions), color='tab:orange')
ax2.axvline(2, color='tab:green', ls='--')

relevant_data = data[data[:, 1] == 450]
dimensions_450 = relevant_data[:, 2]
ax2.hist(dimensions_450, bins=np.arange(0.5, 5.5, 0.25))

plt.savefig(outfile)

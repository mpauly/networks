import matplotlib.pyplot as plt
import numpy as np

dimfile = 'data/dim_roadnet_pa.dat'
pngfile = 'plots/dim_roadnet_pa.png'

print("Reading file {} and writing plot to {}".format(dimfile, pngfile))
# 1D chain plot
data = np.loadtxt(dimfile)
startnode, sigma, dim = data.T
startnodes = np.unique(startnode)

for node in startnodes:
    mask = (node == startnode) & (~np.isnan(dim))
    plt.plot(sigma[mask], dim[mask], alpha=0.2, c='tab:blue')

sigmas = np.unique(sigma)
means = np.zeros(sigmas.shape)
stds = np.zeros(sigmas.shape)

# poor mans implementation of a pivot table
for ind, s in enumerate(sigmas):
    relevant_data = data[data[:, 1] == s]
    dimensions = relevant_data[:, 2]
    means[ind] = np.mean(dimensions)
    stds[ind] = np.std(dimensions)

plt.plot(sigmas, means, c='tab:orange')
plt.fill_between(sigmas, means + stds, means - stds, alpha=0.2, color='tab:orange')

plt.axhline(2, c='tab:green')

plt.xlim([np.min(sigma), np.max(sigma)])
plt.ylim(1, 3)

plt.xlabel('$\\sigma$')
plt.ylabel('$d_{\\rm spec}$')
plt.savefig(pngfile)

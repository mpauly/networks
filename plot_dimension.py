import sys

import matplotlib.pyplot as plt
import numpy as np

if len(sys.argv) != 2:
    raise ValueError('Please provide one dimension file to plot.')

dimfile = sys.argv[1]
pngfile = dimfile.replace('data/', 'plots/').replace('.dat', '.png')

print("Reading file {} and writing plot to {}".format(dimfile, pngfile))
# 1D chain plot
startnode, sigma, dim = np.loadtxt(dimfile, unpack=True)
startnodes = np.unique(startnode)

for node in startnodes:
    mask = (node == startnode) & (~np.isnan(dim))
    plt.plot(sigma[mask], dim[mask])

plt.xlabel('$\\sigma$')
plt.ylabel('$d_{\\rm spec}$')
plt.savefig(pngfile)

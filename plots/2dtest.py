import matplotlib.pyplot as plt
import numpy as np

outfile = "plots/out/2dtest.pdf"

for j in range(9):
    i = 9 - j
    dimfile = "data/dimension/2dtest/d0{}.dat".format(i)
    print("Reading file {} ".format(dimfile))
    # 1D chain plot
    startnode, sigma, dim = np.loadtxt(dimfile, unpack=True)
    startnodes = np.unique(startnode)

    for node in startnodes:
        mask = (node == startnode) & (~np.isnan(dim))
    plt.plot(sigma[mask], dim[mask], label="$\\delta=0.{}$".format(i))

plt.legend()
plt.xlabel("$\\sigma$")
plt.ylabel("$d_{\\rm spec}$")

fig = plt.gcf()
fig.set_size_inches(5.52, 3.41)
fig.savefig(outfile, bbox_inches="tight")

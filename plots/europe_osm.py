import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition

matplotlib.rcParams["mathtext.fontset"] = "stix"
matplotlib.rcParams["font.family"] = "STIXGeneral"

dimfile = "data/dimension/europe_osm_reduced.dat"
outfile = "plots/out/europe_osm.pdf"

fig, ax1 = plt.subplots()

print("Reading file {} and writing plot to {}".format(dimfile, outfile))
# 1D chain plot
data = np.loadtxt(dimfile)
startnode, sigma, dim = data.T
startnodes = np.unique(startnode)

for node in startnodes:
    mask = (node == startnode) & (~np.isnan(dim))
    ax1.plot(sigma[mask], dim[mask], alpha=0.2, c="tab:blue", zorder=0.5)

sigmas = np.unique(sigma)
means = np.zeros(sigmas.shape)
stds = np.zeros(sigmas.shape)

# poor mans implementation of a pivot table
for ind, s in enumerate(sigmas):
    relevant_data = data[data[:, 1] == s]
    dimensions = relevant_data[:, 2]
    means[ind] = np.mean(dimensions)
    stds[ind] = np.std(dimensions, ddof=1)

ax1.plot(sigmas, means, c="tab:orange")
plt.fill_between(sigmas, means + stds, means - stds, alpha=0.2, color="tab:orange")

plt.axhline(2, c="tab:green", ls="--")

plt.xlim([np.min(sigma), np.max(sigma)])
plt.ylim(0, 3.5)

plt.xlabel("$\\sigma$")
plt.ylabel("$d_{\\rm spec}$")

ax2 = plt.axes([0, 0, 1, 1])
# Manually set the position and relative size of the inset axes within ax1
ip = InsetPosition(ax1, [0.65, 0.65, 0.3, 0.3])
sigma_ref = 2500
ax2.set_axes_locator(ip)
ax2.set_xlabel("$d_{{\\rm spec}}(\\sigma = {})$".format(sigma_ref))
ax2.set_ylabel("Freqency")

relevant_data = data[data[:, 1] == sigma_ref]
dimensions_ref = relevant_data[:, 2]
ax2.hist(dimensions_ref, bins=np.arange(0.5, 5.5, 0.25))
ax2.axvline(np.mean(dimensions_ref), color="tab:orange")
ax2.axvline(2, color="tab:green", ls="--")

fig = plt.gcf()
fig.set_size_inches(5.52, 3.41)
fig.savefig(outfile, bbox_inches="tight")

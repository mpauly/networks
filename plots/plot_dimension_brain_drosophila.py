import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition
from scipy.signal import argrelextrema

matplotlib.rcParams["mathtext.fontset"] = "stix"
matplotlib.rcParams["font.family"] = "STIXGeneral"

dimfile = "data/dim_fly-drosophila-large.dat"
outfile = "plots/out/dim_drosophila_large.pdf"

target_dimension = 3
hist_scale = 30

print("Reading file {} and writing plot to {}".format(dimfile, outfile))
# 1D chain plot
data = np.loadtxt(dimfile)
startnode, sigma, dim = data.T
startnodes = np.unique(startnode)

local_min_walkers = []
other_walks = []
late_max_walkers = []

# Sort data into types of trajectories
for node in startnodes:
    mask = (node == startnode) & (~np.isnan(dim))
    local_mins = argrelextrema(dim[mask], np.less)[0]
    global_max = np.argmax(dim[mask])

    walk = {"start_node": node, "sigma": sigma[mask], "dim": dim[mask]}

    has_local_min = local_mins.size > 0 and 10 < local_mins[0] < 50
    if has_local_min and local_mins[0] > global_max:
        local_min_walkers.append(walk)
    elif has_local_min:
        late_max_walkers.append(walk)
    else:
        other_walks.append(walk)

print("Total walker count:")
for walk_set, title in zip(
    [other_walks, late_max_walkers, local_min_walkers],
    [
        "Other walks",
        "Walks with a late maximum",
        "Walks without local minimum and early maximum",
    ],
):
    print("  Walk {}  -  {} walkers".format(title, len(walk_set)))
    print("   Example start node: {}".format(walk_set[-1]["start_node"]))

# plot different types
fig, ax1 = plt.subplots()

for walk_set, color in zip(
    [other_walks, late_max_walkers, local_min_walkers],
    ["tab:green", "tab:orange", "tab:blue"],
):
    for walker in walk_set:
        ax1.plot(walker["sigma"], walker["dim"], alpha=0.2, c=color)

for walk_set, color in zip(
    [other_walks, late_max_walkers, local_min_walkers],
    ["tab:brown", "tab:purple", "tab:red"],
):
    sigmas = np.mean(np.array([w["sigma"] for w in walk_set]), axis=0)
    means = np.mean(np.array([w["dim"] for w in walk_set]), axis=0)
    ax1.plot(sigmas, means, c=color)

plt.axhline(target_dimension, c="tab:gray", ls="--")

plt.xlim([np.min(sigma), np.max(sigma)])
plt.ylim(0.0, 15.0)

plt.xlabel("$\\sigma$")
plt.ylabel("$d_{\\rm spec}$")

ax2 = plt.axes([0, 0, 1, 1])
# Manually set the position and relative size of the inset axes within ax1
ip = InsetPosition(ax1, [0.65, 0.65, 0.3, 0.3])
ax2.set_axes_locator(ip)
ax2.set_xlabel("$d_{{\\rm spec}}(\\sigma = {})$".format(hist_scale))
ax2.set_ylabel("$n$")

all_dims = []

for walk_set in [other_walks, late_max_walkers, local_min_walkers]:
    ref_dimensions = np.array(
        [walker["dim"][walker["sigma"] == hist_scale][0] for walker in walk_set]
    ).flatten()
    all_dims.append(ref_dimensions)

ax2.hist(all_dims, color=["tab:green", "tab:orange", "tab:blue"])
ax2.axvline(target_dimension, color="tab:gray", ls="--")
ax2.set_ylim([0, 100])

fig = plt.gcf()
fig.set_size_inches(5.52, 3.41)
fig.savefig(outfile, bbox_inches="tight")

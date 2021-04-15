# %%
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition
from scipy.signal import argrelextrema
from trajectory_classification import Trajectories

matplotlib.rcParams["mathtext.fontset"] = "stix"
matplotlib.rcParams["font.family"] = "STIXGeneral"


# %%
dimfile = "../data/dimension/fly-drosophila-large.dat"
outfile = "../plots/out/drosophila_large.pdf"

target_dimension = 3
hist_scale = 30

# %%
# 1D chain plot
data = pd.read_table(dimfile, comment="#", names=["start_node", "sigma", "dim"])
# %%
startnodes = pd.unique(data["start_node"])
maxsigma = data["sigma"].max()
data["min"] = data.iloc[argrelextrema(data.dim.values, np.less_equal)[0]]["dim"]
data["max"] = data.iloc[argrelextrema(data.dim.values, np.greater_equal)[0]]["dim"]
data["tratype"] = 0
# %%
data.loc[data["sigma"] == 1, "min"] = np.nan
data.loc[data["sigma"] == 1, "max"] = np.nan
data.loc[data["sigma"] == maxsigma, "min"] = np.nan
data.loc[data["sigma"] == maxsigma, "max"] = np.nan

# %% [markdown]
# Here we are classifying all trajectories into different types
# %%
for sn in startnodes:
    traj = data[data["start_node"] == sn]
    maxima = traj[traj["max"].notnull()]
    data.loc[data["start_node"] == sn, "tratype"] = Trajectories.classify(maxima)
# %%
data["tratype"].value_counts()
# %%
fig, ax1 = plt.subplots()

for walk_type, col in Trajectories.iter():
    data_plot = np.array(
        list(
            data[data["tratype"] == walk_type]
            .groupby("start_node")
            .apply(pd.DataFrame.to_numpy)
        )
    )
    data_plot = data_plot[:, :, 1:3]
    line_segments = LineCollection(data_plot, alpha=0.2, color=col["color"])
    ax1.add_collection(line_segments)

    mean_per_sigma = data[data["tratype"] == walk_type].groupby("sigma").mean()
    ax1.plot(mean_per_sigma.index, mean_per_sigma["dim"], c=col["mean_color"])

plt.axhline(3, c="tab:gray", ls="--")
plt.ylim(0.0, 15.0)

plt.xlabel("$\\sigma$")
plt.ylabel("$d_{\\rm spec}$")

ax2 = plt.axes([0, 0, 1, 1])
# Manually set the position and relative size of the inset axes within ax1
ip = InsetPosition(ax1, [0.65, 0.65, 0.3, 0.3])
ax2.set_axes_locator(ip)
ax2.set_xlabel("$d_{{\\rm spec}}(\\sigma = {})$".format(hist_scale))
ax2.set_ylabel("$n$")

all_dims = [
    data[(data["tratype"] == walk_type) & (data["sigma"] == 30)]["dim"].values
    for walk_type, c in Trajectories.iter()
]
ax2.hist(all_dims, color=[c["color"] for walk_type, c in Trajectories.iter()])

ax2.axvline(target_dimension, color="tab:gray", ls="--")
ax2.set_ylim([0, 100])

fig = plt.gcf()

# %%
fig.set_size_inches(5.52, 3.41)
fig.savefig(outfile, bbox_inches="tight")

# %%

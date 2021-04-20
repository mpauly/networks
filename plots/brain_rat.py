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
# %% [markdown]
# ## Voxel rat brain
# Spectral dimension for a brain graph of a rat.
# Graph is obtained by taking a correlation matrix and introducing
# an arbitrary cutoff.

# %%
dimfile = "../data/dimension/rat_voxel_brain.dat"
outfile = "../plots/out/rat_voxel_brain.pdf"
data = pd.read_table(dimfile, comment="#", names=["start_node", "sigma", "dim"])

data_plot = np.array(list(data.groupby("start_node").apply(pd.DataFrame.to_numpy)))
data_plot = data_plot[:, :, 1:]
# %%
fig, ax1 = plt.subplots()
ax1.set_xlim(data["sigma"].min(), data["sigma"].max())
ax1.set_xlim(0, 110)
line_segments = LineCollection(data_plot, alpha=0.2)
ax1.add_collection(line_segments)
ax1.axhline(3, c="tab:green", ls="--")
ax1.set_xlabel("$\\sigma$")
ax1.set_ylabel("$d_{\\rm spec}$")
# %% [markdown]
# Next we try to apply our classification of trajectories again.
# %%
startnodes = pd.unique(data["start_node"])
maxsigma = data["sigma"].max()
data["min"] = data.iloc[argrelextrema(data.dim.values, np.less_equal, order=3)[0]][
    "dim"
]
data["max"] = data.iloc[argrelextrema(data.dim.values, np.greater_equal, order=3)[0]][
    "dim"
]
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
fig, ax1 = plt.subplots()

for walk_type, col in Trajectories.iter():
    data_plot = np.array(
        list(
            data[data["tratype"] == walk_type]
            .groupby("start_node")
            .apply(pd.DataFrame.to_numpy)
        )
    )
    if data_plot.shape[0] > 0:
        data_plot = data_plot[:, :, 1:3]
        line_segments = LineCollection(data_plot, alpha=0.2, color=col["color"])
        ax1.add_collection(line_segments)

        mean_per_sigma = data[data["tratype"] == walk_type].groupby("sigma").mean()
        ax1.plot(mean_per_sigma.index, mean_per_sigma["dim"], c=col["mean_color"])

plt.axhline(3, c="tab:gray", ls="--")

plt.xlim(0, 110)

plt.xlabel("$\\sigma$")
plt.ylabel("$d_{\\rm spec}$")

# %%
fig.set_size_inches(5.52, 3.41)
fig.savefig(outfile, bbox_inches="tight")
# %%

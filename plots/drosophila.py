# %% [markdown]
# # Setup
# Below we study the *weighted* drosophila graph.

# %%
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition
from scipy.signal import argrelextrema


# %%
matplotlib.rcParams["mathtext.fontset"] = "stix"
matplotlib.rcParams["font.family"] = "STIXGeneral"


# %%
dimfile = "../data/dimension/fly-drosophila-weighted.dat"
outfile = "../plots/out/fly-drosophila-weighted.pdf"


# %%
data = pd.read_table(dimfile, comment="#", names=["start_node", "sigma", "dim"])
startnodes = pd.unique(data["start_node"])
maxsigma = data["sigma"].max()


# %%
startnodes.shape

# %% [markdown]
# ## Standard Plot

# %%
data_plot = np.array(list(data.groupby("start_node").apply(pd.DataFrame.to_numpy)))
data_plot = data_plot[:, :, 1:3]

fig, ax1 = plt.subplots()
line_segments = LineCollection(data_plot, alpha=0.2)
ax1.add_collection(line_segments)

mean_per_sigma = data.groupby("sigma").mean()
std_dev_per_sigma = data.groupby("sigma").agg(np.std, ddof=1)
print(
    "Maximum is at {} and has dimension {}".format(
        mean_per_sigma["dim"].idxmax(), mean_per_sigma["dim"].max()
    )
)

mean_per_sigma["dim"].plot(ax=ax1, color="tab:orange")
plt.fill_between(
    mean_per_sigma.index,
    mean_per_sigma["dim"] + std_dev_per_sigma["dim"],
    mean_per_sigma["dim"] - std_dev_per_sigma["dim"],
    alpha=0.2,
    color="tab:orange",
)

ax1.set_xlim(data["sigma"].min(), data["sigma"].max())
ax1.set_ylim(0, 12)
# plt.axhline(2, c="tab:green", ls="--")
plt.xlabel("$\\sigma$")
plt.ylabel("$d_{\\rm spec}$")

# %% [markdown]
# ## Splitting into various types of trajectories

# %%
data["min"] = data.iloc[argrelextrema(data.dim.values, np.less_equal, order=3)[0]][
    "dim"
]
data["max"] = data.iloc[argrelextrema(data.dim.values, np.greater_equal, order=3)[0]][
    "dim"
]
data["tratype"] = 0
data.loc[data["sigma"] == 1, "min"] = np.nan
data.loc[data["sigma"] == 1, "max"] = np.nan
data.loc[data["sigma"] == maxsigma, "min"] = np.nan
data.loc[data["sigma"] == maxsigma, "max"] = np.nan


# %%
data["tratype"] = 0
for sn in startnodes:
    traj = data[data["start_node"] == sn]
    maxima = traj[traj["max"].notnull()]
    # trajectory has less then two maxima
    if maxima.shape[0] < 2:
        data.loc[data["start_node"] == sn, "tratype"] = 3
        continue
    # first maximum is higher than second maximum
    if maxima.iloc[0]["dim"] > maxima.iloc[1]["dim"]:
        data.loc[data["start_node"] == sn, "tratype"] = 1
    # second is higher than first
    if maxima.iloc[0]["dim"] < maxima.iloc[1]["dim"]:
        data.loc[data["start_node"] == sn, "tratype"] = 2


# %%
target_dimension = 3
fig, ax = plt.subplots(3, figsize=(5, 12))
for ttype in range(3):
    mean_per_sigma = data[data["tratype"] == ttype + 1].groupby("sigma").mean()
    data_plot = np.array(
        list(
            data[data["tratype"] == ttype + 1]
            .groupby("start_node")
            .apply(pd.DataFrame.to_numpy)
        )
    )
    data_plot = data_plot[:, :, 1:3]
    line_segments = LineCollection(data_plot, alpha=0.2)
    ax[ttype].add_collection(line_segments)
    ax[ttype].plot(mean_per_sigma.index, mean_per_sigma["dim"], c="tab:red")
    ax[ttype].axhline(target_dimension, color="tab:gray", ls="--")


# %%

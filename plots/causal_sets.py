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
dimfile = "../data/dimension/causal_set_l30.dat"

# %%
# 1D chain plot
data = pd.read_table(dimfile, comment="#", names=["start_node", "sigma", "dim"])
startnodes = pd.unique(data["start_node"])
maxsigma = data["sigma"].max()
# %%
data_plot = np.array(list(data.groupby("start_node").apply(pd.DataFrame.to_numpy)))
data_plot = data_plot[:, :, 1:3]

fig, ax1 = plt.subplots()
ax1.set_xlim(data["sigma"].min(), data["sigma"].max())
ax1.set_ylim(0, 15)
line_segments = LineCollection(data_plot, alpha=0.2)
ax1.add_collection(line_segments)

mean_per_sigma = data.groupby("sigma").mean()
std_dev_per_sigma = data.groupby("sigma").agg(np.std, ddof=1)
print(
    "Maximum is at {} and has dimension {}".format(
        mean_per_sigma["dim"].idxmax(), mean_per_sigma["dim"].max()
    )
)
ref_sigma = mean_per_sigma["dim"].idxmax()

mean_per_sigma["dim"].plot(ax=ax1, color="tab:orange")
plt.fill_between(
    mean_per_sigma.index,
    mean_per_sigma["dim"] + std_dev_per_sigma["dim"],
    mean_per_sigma["dim"] - std_dev_per_sigma["dim"],
    alpha=0.2,
    color="tab:orange",
)

# plt.axhline(2, c="tab:green", ls="--")
plt.xlabel("$\\sigma$")
plt.ylabel("$d_{\\rm spec}$")


# %%
fig.set_size_inches(5.52, 3.41)
fig.savefig("../plots/out/causal_set_l30.pdf", bbox_inches="tight")

# %%

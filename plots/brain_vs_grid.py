# %%
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition

matplotlib.rcParams["mathtext.fontset"] = "stix"
matplotlib.rcParams["font.family"] = "STIXGeneral"
# %% [markdown]
# ## Human Brain vs Grid
# Below we compare the spectral dimension of a 3D grid
# to the brain split into roughly 300 regions.

# %%
dimfile = "../data/dimension/human_brain_regions_300.dat"
outfile = "../plots/out/human_brain_regions_300.pdf"
data = pd.read_table(dimfile, comment="#", names=["start_node", "sigma", "dim"])

data_plot = np.array(list(data.groupby("start_node").apply(pd.DataFrame.to_numpy)))
data_plot = data_plot[:, :, 1:]
# %%
fig, ax1 = plt.subplots()
ax1.set_xlim(data["sigma"].min(), data["sigma"].max())
ax1.set_ylim(0, 10)
ax1.set_xlim(0, 30)
line_segments = LineCollection(data_plot, alpha=0.2)
ax1.add_collection(line_segments)
ax1.axhline(3, c="tab:green", ls="--")
ax1.set_xlabel("$\\sigma$")
ax1.set_ylabel("$d_{\\rm spec}$")
#%%
fig.set_size_inches(5.52, 3.41)
fig.savefig(outfile, bbox_inches="tight")
# %% [markdown]
# The brain is extremly well connected, the spectral dimension
# decays more or less instantly.
# %%
dimfile = "../data/dimension/3d_lattice_7.dat"
outfile = "../plots/out/3d_lattice_7.pdf"
data = pd.read_table(dimfile, comment="#", names=["start_node", "sigma", "dim"])

data_plot = np.array(list(data.groupby("start_node").apply(pd.DataFrame.to_numpy)))
data_plot = data_plot[:, :, 1:]
# %%
plt.clf()
fig, ax1 = plt.subplots()
ax1.set_xlim(data["sigma"].min(), data["sigma"].max())
ax1.set_ylim(0, 4)
line_segments = LineCollection(data_plot, alpha=0.2)
ax1.add_collection(line_segments)

ax1.axhline(3, c="tab:green", ls="--")
ax1.set_xlabel("$\\sigma$")
ax1.set_ylabel("$d_{\\rm spec}$")
#%%
fig.set_size_inches(5.52, 3.41)
fig.savefig(outfile, bbox_inches="tight")
# %% [markdown]
# Various of the walkers actually experience a dimension smaller than three.
# Presumably due to the finite extend.
# %%
dimfile = "../data/dimension/3d_lattice_100.dat"
data = pd.read_table(dimfile, comment="#", names=["start_node", "sigma", "dim"])

data_plot = np.array(list(data.groupby("start_node").apply(pd.DataFrame.to_numpy)))
data_plot = data_plot[:, :, 1:]
# %%
plt.clf()
fig, ax1 = plt.subplots()
ax1.set_xlim(data["sigma"].min(), data["sigma"].max())
ax1.set_ylim(0, 10)
line_segments = LineCollection(data_plot, alpha=0.2)
ax1.add_collection(line_segments)

ax1.axhline(3, c="tab:green", ls="--")
ax1.set_xlabel("$\\sigma$")
ax1.set_ylabel("$d_{\\rm spec}$")
# %%

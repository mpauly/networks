# %% [markdown]
# ## Watts Strogatz
# %%
from itertools import cycle

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition

matplotlib.rcParams["mathtext.fontset"] = "stix"
matplotlib.rcParams["font.family"] = "STIXGeneral"
# %%
dimfile = "../data/dimension/watts_strogatz.dat"
data = pd.read_table(dimfile)

# %% [markdown]
# Different diffusion constants do not matter for the general finding

# %%
ref_degree = 4
format_din = (11.04, 3.41)
fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(15, 4))
for ind, diffconst in enumerate([0.25, 0.5, 0.75]):
    for name, group in data[
        (data["diffusion_const"] == diffconst) & (data["degree"] == ref_degree)
    ].groupby("rewiring_prob"):
        ax[ind].plot(group["sigma"], group["dim"], label="$\\beta={}$".format(name))
    ax[ind].legend()
    ax[ind].set_xlim(0, 200)
    ax[ind].set_xlabel("$\\sigma$")
    ax[ind].set_ylabel("$d_{\\rm spec}$")
    ax[ind].set_title("$\\delta={}$".format(diffconst))

# %% [markdown]
# Due to the effectively 2d structure one obtains the intermediate dip

# %%
fig.clf()
ax = fig.gca()
for name, group in data[
    (data["diffusion_const"] == 0.5) & (data["degree"] == ref_degree)
].groupby("rewiring_prob"):
    ax.plot(group["sigma"], group["dim"], label="$\\beta={}$".format(name))
ax.legend()
ax.set_xlim(0, 200)
ax.set_xlabel("$\\sigma$")
ax.set_ylabel("$d_{\\rm spec}$")
ax.set_title("$\\delta={}$".format(0.5))
fig


# %%
lines = ["-", "--", "-.", ":"]
ref_delta = 0.5
format_din = (8.04, 2.41)
fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(15, 4))
for ind, degree in enumerate([2, 4, 6]):
    linecycler = cycle(lines)
    for name, group in data[
        (data["diffusion_const"] == ref_delta) & (data["degree"] == degree)
    ].groupby("rewiring_prob"):
        ax[ind].plot(
            group["sigma"],
            group["dim"],
            label="$\\beta={}$".format(name),
            ls=next(linecycler),
        )
    ax[ind].legend()
    ax[ind].set_xlim(0, 200)
    ax[ind].set_ylim(0, 7)
    ax[ind].set_xlabel("$\\sigma$")
    ax[ind].set_ylabel("$d_{\\rm spec}$")
    ax[ind].set_title("$\\bar{k}=" + str(degree) + "$")

# %%
fig.set_size_inches(11.04, 3.41)
fig.savefig("../plots/out/watts_strogatz.pdf", bbox_inches="tight")
# %% [markdown]
# ## 2D generalization
# Here $\beta=0.001$, i.e. we are in the case of small deviations from a grid
# %%
dimfile = "../data/dimension/watts_strogatz_2d_100.dat"
data = pd.read_table(dimfile, comment="#", names=["start_node", "sigma", "dim"])
data_plot = np.array(list(data.groupby("start_node").apply(pd.DataFrame.to_numpy)))
data_plot = data_plot[:, :, 1:]
#%%
mean_per_sigma = data.groupby("sigma").mean()
std_dev_per_sigma = data.groupby("sigma").agg(np.std, ddof=1)
print(
    "Maximum is at {} and has dimension {}".format(
        mean_per_sigma["dim"].idxmax(), mean_per_sigma["dim"].max()
    )
)
ref_sigma = mean_per_sigma["dim"].idxmax()
# %%
plt.clf()
fig, ax1 = plt.subplots()
ax1.set_xlim(data["sigma"].min(), data["sigma"].max())
line_segments = LineCollection(data_plot, alpha=0.2)
ax1.add_collection(line_segments)

mean_per_sigma["dim"].plot(ax=ax1, color="tab:orange")
ax1.fill_between(
    mean_per_sigma.index,
    mean_per_sigma["dim"] + std_dev_per_sigma["dim"],
    mean_per_sigma["dim"] - std_dev_per_sigma["dim"],
    alpha=0.2,
    color="tab:orange",
)

ax1.axhline(2, c="tab:green", ls="--")
ax1.set_xlabel("$\\sigma$")
ax1.set_ylabel("$d_{\\rm spec}$")

ax2 = plt.axes([0, 0, 1, 1])
ip = InsetPosition(ax1, [0.65, 0.65, 0.3, 0.3])
ax2.set_axes_locator(ip)

line_segments = LineCollection(data_plot, alpha=0.2)
ax2.add_collection(line_segments)
mean_per_sigma["dim"].plot(ax=ax2, color="tab:orange")
ax2.fill_between(
    mean_per_sigma.index,
    mean_per_sigma["dim"] + std_dev_per_sigma["dim"],
    mean_per_sigma["dim"] - std_dev_per_sigma["dim"],
    alpha=0.2,
    color="tab:orange",
)
ax2.set_xlabel("$\\sigma$")
ax2.set_ylabel("$d_{\\rm spec}$")
ax2.set_xlim(0, 25)
# %%
fig.set_size_inches(5.52, 3.41)
fig.savefig("../plots/out/watts_strogatz_2d.pdf", bbox_inches="tight")
# %% [markdown]
# ## 3D generalization
# %%
dimfile = "../data/dimension/watts_strogatz_3d_100_0.001000.dat"
data = pd.read_table(dimfile, comment="#", names=["start_node", "sigma", "dim"])
data_plot = np.array(list(data.groupby("start_node").apply(pd.DataFrame.to_numpy)))
data_plot = data_plot[:, :, 1:]
# %%
mean_per_sigma = data.groupby("sigma").mean()
std_dev_per_sigma = data.groupby("sigma").agg(np.std, ddof=1)
print(
    "Maximum is at {} and has dimension {}".format(
        mean_per_sigma["dim"].idxmax(), mean_per_sigma["dim"].max()
    )
)
ref_sigma = mean_per_sigma["dim"].idxmax()
#%%
plt.clf()
fig, ax1 = plt.subplots()
ax1.set_xlim(data["sigma"].min(), data["sigma"].max())
line_segments = LineCollection(data_plot, alpha=0.2)
ax1.add_collection(line_segments)

mean_per_sigma["dim"].plot(ax=ax1, color="tab:orange")
ax1.fill_between(
    mean_per_sigma.index,
    mean_per_sigma["dim"] + std_dev_per_sigma["dim"],
    mean_per_sigma["dim"] - std_dev_per_sigma["dim"],
    alpha=0.2,
    color="tab:orange",
)

ax1.axhline(3, c="tab:green", ls="--")
ax1.set_xlabel("$\\sigma$")
ax1.set_ylabel("$d_{\\rm spec}$")

ax2 = plt.axes([0, 0, 1, 1])
ip = InsetPosition(ax1, [0.65, 0.65, 0.3, 0.3])
ax2.set_axes_locator(ip)

line_segments = LineCollection(data_plot, alpha=0.2)
ax2.add_collection(line_segments)
mean_per_sigma["dim"].plot(ax=ax2, color="tab:orange")
ax2.fill_between(
    mean_per_sigma.index,
    mean_per_sigma["dim"] + std_dev_per_sigma["dim"],
    mean_per_sigma["dim"] - std_dev_per_sigma["dim"],
    alpha=0.2,
    color="tab:orange",
)

ax2.set_xlabel("$\\sigma$")
ax2.set_ylabel("$d_{\\rm spec}$")
ax2.set_xlim(0, 25)


# %%
fig.set_size_inches(5.52, 3.41)
fig.savefig("../plots/out/watts_strogatz_3d.pdf", bbox_inches="tight")

#%%

plt.clf()
# %%
lines = ["-", "--", "-.", ":"]
linecycler = cycle(lines)

fig = plt.gcf()
for beta in ["0.001", "0.005", "0.010", "0.050"]:
    dimfile = "../data/dimension/watts_strogatz_3d_100_{}000.dat".format(beta)
    data = pd.read_table(dimfile, comment="#", names=["start_node", "sigma", "dim"])
    mean_per_sigma = data.groupby("sigma").mean()
    mean_per_sigma["dim"].plot(label="$\\beta={}$".format(beta), ls=next(linecycler))

plt.legend()
plt.xlabel("$\\sigma$")
plt.ylabel("$d_{\\rm spec}$")

plt.axhline(3.0, c="tab:gray", ls="--")
# %%
fig.set_size_inches(4.5, 2.7)
fig.savefig("../plots/out/watts_strogatz_3d_var_beta.pdf", bbox_inches="tight")

# %%

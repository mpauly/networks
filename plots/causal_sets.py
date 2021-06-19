# %%
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition
from scipy.signal import argrelextrema
from sklearn import linear_model

# %%
matplotlib.rcParams["mathtext.fontset"] = "stix"
matplotlib.rcParams["font.family"] = "STIXGeneral"

# %%
dfs = []
for dim in [2, 3, 4, 5, 6]:
    thisdata = pd.read_table("../data/other/average_path_{}d.tsv".format(dim))
    thisdata["dim"] = dim
    dfs.append(thisdata)

data = pd.concat(dfs)
data["stderr_spl"] = data["stddev_spl"] / np.sqrt(0.1 * data["nodes"])
data["log_nodes"] = np.log(data["nodes"])
# %%

ax = plt.gca()
dim_df = pd.DataFrame(columns=["dim", "intercept", "slope"])

for dim, col in zip(
    [2, 3, 4, 5, 6], ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple"]
):
    thisdata = data[data["dim"] == dim]
    thisdata.plot.scatter(
        "nodes",
        "avg_spl",
        yerr="stderr_spl",
        ax=ax,
        c=col,
        label="$d={}$".format(dim),
        figsize=(5.52, 3.41),
    )
    regr = linear_model.LinearRegression()
    model = regr.fit(thisdata[["log_nodes"]], thisdata[["avg_spl"]])

    prediction = model.predict(thisdata[["log_nodes"]])
    ax.plot(thisdata[["nodes"]], prediction, c=col)

    dim_df = dim_df.append(
        {"dim": dim, "intercept": model.intercept_[0], "slope": model.coef_[0, 0]},
        ignore_index=True,
    )

plt.legend()
plt.xlabel("Nr. of nodes")
plt.ylabel("Mean shortest path length")
plt.semilogx()

# %% [markdown]
# ## Causal Set without regulator
# %%
dimfile = "../data/dimension/causal_set.dat"

# %%
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

plt.axhline(2, c="tab:green", ls="--")

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
fig.savefig("../plots/out/causal_set_noreg.pdf", bbox_inches="tight")


# %% [markdown]
# ## Regulated Causal Sets
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

plt.axhline(2, c="tab:green", ls="--")

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

# %% [markdown]
# # Transitive Percolations
# %%
dimfile = "../data/dimension/transitive_percolations_0.000001.dat"
data = pd.read_table(dimfile, comment="#", names=["start_node", "sigma", "dim"])
startnodes = pd.unique(data["start_node"])
maxsigma = data["sigma"].max()
# %%
data_plot = np.array(list(data.groupby("start_node").apply(pd.DataFrame.to_numpy)))
data_plot = data_plot[:, :, 1:3]

fig, ax1 = plt.subplots()
ax1.set_xlim(data["sigma"].min(), data["sigma"].max())
ax1.set_ylim(0, 8)
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
plt.xlabel("$\\sigma$")
plt.ylabel("$d_{\\rm spec}$")
# %%
fig.set_size_inches(5.52, 3.41)
fig.savefig("../plots/out/transitive_percolations_1.pdf", bbox_inches="tight")
# %%
dimfile = "../data/dimension/transitive_percolations_0.000005.dat"
data = pd.read_table(dimfile, comment="#", names=["start_node", "sigma", "dim"])
startnodes = pd.unique(data["start_node"])
maxsigma = data["sigma"].max()
# %%
data_plot = np.array(list(data.groupby("start_node").apply(pd.DataFrame.to_numpy)))
data_plot = data_plot[:, :, 1:3]

fig, ax1 = plt.subplots()
ax1.set_xlim(data["sigma"].min(), data["sigma"].max())
ax1.set_ylim(0, 11)
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
plt.xlabel("$\\sigma$")
plt.ylabel("$d_{\\rm spec}$")
# %%
fig.set_size_inches(5.52, 3.41)
fig.savefig("../plots/out/transitive_percolations_2.pdf", bbox_inches="tight")
# %%

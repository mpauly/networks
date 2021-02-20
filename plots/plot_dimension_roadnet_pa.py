import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition

matplotlib.rcParams["mathtext.fontset"] = "stix"
matplotlib.rcParams["font.family"] = "STIXGeneral"

dimfile = "data/dim_roadnet_pa.dat"
outfile = "plots/out/dim_roadnet_pa.pdf"
outfile_log = "plots/out/dim_roadnet_pa_log.pdf"
convergence_file = "plots/out/dim_roadnet_pa_convergence.pdf"
convergence_file_2 = "plots/out/dim_roadnet_pa_convergence_2.pdf"

print("Reading file {} and writing plot to {}".format(dimfile, outfile))
# 1D chain plot
data = pd.read_table(dimfile, comment="#", names=["start_node", "sigma", "dim"])
startnodes = pd.unique(data["start_node"])

data_plot = np.array(list(data.groupby("start_node").apply(pd.DataFrame.to_numpy)))
data_plot = data_plot[:, :, 1:]

fig, ax1 = plt.subplots()
ax1.set_xlim(data["sigma"].min(), data["sigma"].max())
ax1.set_ylim(0, 5)
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

plt.axhline(2, c="tab:green", ls="--")
plt.xlabel("$\\sigma$")
plt.ylabel("$d_{\\rm spec}$")

ax2 = plt.axes([0, 0, 1, 1])
# Manually set the position and relative size of the inset axes within ax1
ip = InsetPosition(ax1, [0.65, 0.65, 0.3, 0.3])
ax2.set_axes_locator(ip)
ax2.set_xlabel("$d_{\\rm spec}(\\sigma = " + str(round(ref_sigma)) + ")$")
ax2.set_ylabel("$n$")

relevant_data = data[data["sigma"] == ref_sigma]
relevant_data["dim"].plot(
    kind="hist", ax=ax2, bins=np.arange(0.5, 5.5, 0.25), grid=False
)
ax2.axvline(relevant_data["dim"].mean(), color="tab:orange")
ax2.axvline(2, color="tab:green", ls="--")

fig = plt.gcf()
fig.set_size_inches(5.52, 3.41)
fig.savefig(outfile, bbox_inches="tight")

ax1.semilogx()
fig.savefig(outfile_log, bbox_inches="tight")

plt.clf()

# ------------- convergence plot ----------------------

nr_of_bins = 4

per_bucket = int(data.shape[0] / nr_of_bins)

ax1 = plt.gca()

for quartile in range(nr_of_bins):
    quartile_data = data[quartile * per_bucket : (quartile + 1) * per_bucket]
    mean_per_sigma = quartile_data.groupby("sigma").mean()
    std_dev_per_sigma = quartile_data.groupby("sigma").agg(np.std)

    mean_per_sigma["dim"].plot(ax=ax1, ls="--")
    plt.fill_between(
        mean_per_sigma.index,
        mean_per_sigma["dim"] + std_dev_per_sigma["dim"],
        mean_per_sigma["dim"] - std_dev_per_sigma["dim"],
        alpha=0.1,
        color=ax1.lines[-1].get_color(),
    )


mean_per_sigma = data.groupby("sigma").mean()
std_dev_per_sigma = data.groupby("sigma").agg(np.std)

mean_per_sigma["dim"].plot(ax=ax1)
plt.fill_between(
    mean_per_sigma.index,
    mean_per_sigma["dim"] + std_dev_per_sigma["dim"],
    mean_per_sigma["dim"] - std_dev_per_sigma["dim"],
    alpha=0.1,
    color=ax1.lines[-1].get_color(),
)

plt.axhline(2, c="tab:gray", ls="--")

ax1.set_xlim(data["sigma"].min(), data["sigma"].max())
ax1.set_ylim(1, 3)

plt.xlabel("$\\sigma$")
plt.ylabel("$d_{\\rm spec}$")

fig = plt.gcf()
fig.set_size_inches(5.52, 3.41)
fig.savefig(convergence_file, bbox_inches="tight")

plt.clf()

# ------------- convergence plot 2 ----------------------

sigmas_per_startnode = data[data["start_node"] == 950557].shape[0]

data["diff"] = data.dim - data.groupby("sigma").transform("mean").dim
data["reldiff"] = data["diff"] / data.groupby("sigma").transform("mean").dim
data["reldiff"] = data["reldiff"].apply(np.abs)

rolling_mean_data = data.groupby("sigma")["dim"].expanding().mean().reset_index()
rolling_mean_data["startid"] = (
    rolling_mean_data["level_1"] / sigmas_per_startnode
).apply(np.floor)
del rolling_mean_data["level_1"]
data_pivot = rolling_mean_data.pivot("sigma", "startid")

x = data_pivot.columns.levels[1].values
y = data_pivot.index.values
z = data_pivot.values
xi, yi = np.meshgrid(x, y)
cs = plt.contourf(
    yi, xi, z, alpha=0.7, cmap=plt.cm.jet, levels=np.linspace(1.5, 2.5, 11)
)

plt.xlabel("$\\sigma$")
plt.ylabel("$N$")

fig = plt.gcf()
cbar = fig.colorbar(cs)
cbar.ax.set_ylabel("$\\bar{d}_{\\rm spec}(N)$", rotation=90)
fig.set_size_inches(5.52, 3.41)
fig.savefig(convergence_file_2, bbox_inches="tight")

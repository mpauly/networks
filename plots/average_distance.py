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

radii_file = "../data/radii/internet_caida_as.dat"
rfile = open(radii_file, "r")
# %%
nr_of_lines = 0

list_data = []
for line in rfile.readlines():
    if line[0] == "#" or line[0] == "\n":
        continue
    line = line.split("\t")
    start_node = int(line[0])
    sigma = int(line[1])
    current_radius = len(line) - 2
    data_frame_lines = [
        [start_node, sigma, rad, float(val.strip())] for rad, val in enumerate(line[2:])
    ]
    list_data.extend(data_frame_lines)
    nr_of_lines += 1

data = pd.DataFrame(list_data, columns=["start_node", "sigma", "radius", "prob"])
print("{} lines read".format(nr_of_lines))
# %%
data["weighted_r"] = data["radius"] * data["prob"]
# %%
ax = plt.gca()
data_plot = (
    data.groupby(["start_node", "sigma"])
    .sum("weighted_r")
    .reset_index()
    .apply(pd.DataFrame.to_numpy)
)
data_plot = np.array(list(data_plot.groupby("start_node").apply(pd.DataFrame.to_numpy)))
data_plot = data_plot[:, :, (1, -1)]
line_segments = LineCollection(data_plot, alpha=0.5)
ax.add_collection(line_segments)
plt.xlim(20, 100)
plt.ylim(2, 6)
# %% [markdown]
# ## Small world scaling
# %%
# %%

radii_file = "../data/radii/internet_caida_as_nodes.dat"
rfile = open(radii_file, "r")
# %%
nr_of_lines = 0

list_data = []
for line in rfile.readlines():
    if line[0] == "#" or line[0] == "\n":
        continue
    line = line.split("\t")
    start_node = int(line[0])
    data_frame_lines = [
        [start_node, rad, int(cnt.strip())] for rad, cnt in enumerate(line[2:])
    ]
    list_data.extend(data_frame_lines)
    nr_of_lines += 1

data = pd.DataFrame(list_data, columns=["start_node", "radius", "count"])
# %%
ax = plt.gca()
plot_data = list(data.groupby("start_node").apply(pd.DataFrame.to_numpy))
plot_data = [p[:, 1:3] for p in plot_data]
line_segments = LineCollection(plot_data, alpha=0.5)
ax.add_collection(line_segments)
plt.xlim(0, 9)
plt.yscale("log")
plt.ylim(1, 50000)
# %%

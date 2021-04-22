# %%
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition

matplotlib.rcParams["mathtext.fontset"] = "stix"
matplotlib.rcParams["font.family"] = "STIXGeneral"
# %%
data_nodes = np.loadtxt("../data/other/nodes_cs_minkowski.csv")
data_edges = np.loadtxt("../data/other/edges_cs_minkowski.csv")

nr_of_nodes = data_nodes.shape[0]
data_edges = data_edges[:, 0:4]

# Rescale to get nice axis labels
data_nodes = np.sqrt(nr_of_nodes) * data_nodes
data_edges = np.sqrt(nr_of_nodes) * data_edges

# %%
# the y direction is the time direction, i.e. the first entry
plt.scatter(data_nodes[:, 2], data_nodes[:, 1], s=2, alpha=0.5, color="tab:blue")
# %%
# for label, t, x in data_nodes.tolist():
#    plt.annotate(str(int(label)), (x + 0.025, t - 0.01))

ax = plt.gca()
# reshape to t,x pairs
data_edges = data_edges.reshape(data_edges.shape[0], 2, 2)
# reverse to x,t pairs for plotting
data_edges = data_edges[:, :, ::-1]
line_segments = LineCollection(data_edges, alpha=0.3, linewidths=0.5, color="tab:blue")
ax.add_collection(line_segments)
# %%
ref_node = 250
coords = data_nodes[ref_node, 1:3]
coords = coords[::-1]
ingoing_edges = data_edges[np.all(np.isclose(data_edges[:, 1], coords), axis=1)]
outgoing_edges = data_edges[np.all(np.isclose(data_edges[:, 0], coords), axis=1)]

line_segments = LineCollection(
    ingoing_edges, alpha=0.75, linewidths=1, color="tab:green"
)
ax.add_collection(line_segments)
line_segments = LineCollection(
    outgoing_edges, alpha=0.75, linewidths=1, color="tab:red"
)
ax.add_collection(line_segments)

plt.scatter(
    data_nodes[ref_node, 2], data_nodes[ref_node, 1], s=6, zorder=2.5, color="black"
)
# %%
plt.xlim(-0.1, np.sqrt(nr_of_nodes) + 0.1)
plt.ylim(-0.1, np.sqrt(nr_of_nodes) + 0.1)
plt.xlabel("$x$")
plt.ylabel("$t$")
plt.tight_layout()
fig = plt.gcf()
fig.set_size_inches(3.41, 3.41)
fig.savefig("../plots/out/causal_set_example.pdf", bbox_inches="tight")

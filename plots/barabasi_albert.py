# %% [markdown]
# ## Barabassi Albert
# Small network, $10^3$ nodes
# %%
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

matplotlib.rcParams["mathtext.fontset"] = "stix"
matplotlib.rcParams["font.family"] = "STIXGeneral"
# %%
dimfile = "../data/dimension/barabasi_albert.dat"
data = pd.read_table(dimfile)

# %%
format_din = (11.04, 3.41)
fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(15, 4))
for ind, diffconst in enumerate([0.25, 0.5, 0.75]):
    for name, group in data[(data["diffusion_const"] == diffconst)].groupby("degree"):
        ax[ind].plot(group["sigma"], group["dim"], label="$D={}$".format(name))
    ax[ind].legend()
    ax[ind].set_xlim(0, 200)
    ax[ind].set_xlabel("$\\sigma$")
    ax[ind].set_ylabel("$d_{\\rm spec}$")
    ax[ind].set_title("$\\delta={}$".format(diffconst))

# %%
fig.set_size_inches(11.04, 3.41)
fig.savefig("../plots/out/barabasi_albert_small.pdf", bbox_inches="tight")

# %% [markdown]
# ## Barabassi Albert - large version
# $N=10^5$ nodes
# %%
dimfile = "../data/dimension/barabasi_albert_large.dat"
data_large = pd.read_table(dimfile)

# %%
format_din = (11.04, 3.41)
fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(15, 4))
for ind, diffconst in enumerate([0.25, 0.5, 0.75]):
    for name, group in data_large[(data_large["diffusion_const"] == diffconst)].groupby(
        "degree"
    ):
        ax[ind].plot(group["sigma"], group["dim"], label="$D={}$".format(name))
    ax[ind].legend()
    ax[ind].set_xlim(0, 200)
    ax[ind].set_xlabel("$\\sigma$")
    ax[ind].set_ylabel("$d_{\\rm spec}$")
    ax[ind].set_title("$\\delta={}$, large".format(diffconst))

# %%
fig.set_size_inches(11.04, 3.41)
fig.savefig("../plots/out/barabasi_albert_large.pdf", bbox_inches="tight")

# %% [markdown]
# ## Barabassi Albert - Combination
# %%
plt.clf()

format_din = (11.04, 3.41)
fig, ax = plt.subplots(nrows=1, ncols=2, figsize=format_din)
diffconst = 0.5
sizes = [3, 5]
for ind, d in enumerate([data, data_large]):
    for name, group in d[(d["diffusion_const"] == diffconst)].groupby("degree"):
        ax[ind].plot(group["sigma"], group["dim"], label="$D={}$".format(name))
    ax[ind].legend()
    ax[ind].set_xlim(0, 200)
    ax[ind].set_xlabel("$\\sigma$")
    ax[ind].set_ylabel("$d_{\\rm spec}$")
    ax[ind].set_title("$N=10^{}$".format(sizes[ind]))
# %%
fig.set_size_inches(11.04, 3.41)
fig.savefig("../plots/out/barabasi_albert_various_size.pdf", bbox_inches="tight")

# %%

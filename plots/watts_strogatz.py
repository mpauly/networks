# %% [markdown]
# ## Setup
# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# %%
dimfile = "../data/dimension/watts_strogatz.dat"


# %%
data = pd.read_table(dimfile)

# %% [markdown]
# Different diffusion constants do not matter for the general finding

# %%
format_din = (11.04, 3.41)
fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(15, 4))
for ind, diffconst in enumerate([0.25, 0.5, 0.75]):
    for name, group in data[data["diffusion_const"] == diffconst].groupby(
        "rewiring_prob"
    ):
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
for name, group in data[data["diffusion_const"] == 0.5].groupby("rewiring_prob"):
    ax.plot(group["sigma"], group["dim"], label="$\\beta={}$".format(name))
ax.legend()
ax.set_xlim(0, 200)
ax.set_xlabel("$\\sigma$")
ax.set_ylabel("$d_{\\rm spec}$")
ax.set_title("$\\delta={}$".format(0.5))
fig


# %%

from itertools import cycle

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition

matplotlib.rcParams["mathtext.fontset"] = "stix"
matplotlib.rcParams["font.family"] = "STIXGeneral"

dim_dir = "data/dimension/"

dimfile = dim_dir + "ring_random.dat"

reference_delta = 0.5
reference_conn = 200

lines = ["-", "--", "-.", ":"]
linecycler = cycle(lines)

print("Reading file {}".format(dimfile))

nr_conn, diffusion_const, sigma, dim = np.loadtxt(dimfile, unpack=True)
nr_conn_values = np.unique(nr_conn)
diffusion_const_values = np.unique(diffusion_const)

# first plot
for nc in nr_conn_values:
    mask = (nc == nr_conn) & np.equal(diffusion_const, reference_delta)
    plt.plot(
        sigma[mask],
        dim[mask],
        label="${:d}, \\bar{{k}} = {:1.1f}$".format(
            round(nc), 2 * (round(nc) + 100) / 100
        ),
        ls=next(linecycler),
    )

plt.legend()
# plt.title("$\\delta = {}$, varying number of connections".format(reference_delta))
plt.xlabel("$\\sigma$")
plt.ylabel("$d_{\\rm spec}$")
plt.xlim(0, 200)

fig = plt.gcf()
fig.set_size_inches(5.52, 3.41)
fig.savefig("plots/out/ring_1.pdf", bbox_inches="tight")

plt.clf()

# second plot
for dc in diffusion_const_values:
    mask = (dc == diffusion_const) & np.equal(nr_conn, reference_conn)
    plt.plot(sigma[mask], dim[mask], label="$\\delta={:.2f}$".format(dc))

plt.legend()
plt.title("${}$ connections".format(reference_conn))
plt.xlabel("$\\sigma$")
plt.ylabel("$d_{\\rm spec}$")
plt.xlim(0, 200)

fig = plt.gcf()
fig.set_size_inches(5.52, 3.41)
fig.savefig("plots/out/ring_2.pdf", bbox_inches="tight")

plt.clf()

# third plot

dimfile = dim_dir + "ring.dat"
print("Reading file {}".format(dimfile))
diffusion_const, sigma, dim = np.loadtxt(dimfile, unpack=True)
diffusion_const_values = np.unique(diffusion_const)

chain_length = 100

for dc in diffusion_const_values:
    mask = dc == diffusion_const
    these_sigmas = sigma[mask]
    these_dim = dim[mask]
    line = plt.plot(
        these_sigmas[:-2], these_dim[:-2], label="$\\delta={:.2f}$".format(dc)
    )

    eigenvalue1 = dc * (1.0 - np.cos(2.0 * np.pi * 1.0 / chain_length))
    theoretical_expectation = (
        4.0 * these_sigmas * eigenvalue1 * np.exp(-eigenvalue1 * these_sigmas)
    )
    plt.plot(
        these_sigmas,
        theoretical_expectation,
        ls="--",
        color=line[0].get_color(),
        alpha=0.5,
    )

plt.legend()
plt.title("${}$ connections".format(0))
plt.xlabel("$\\sigma$")
plt.ylabel("$d_{\\rm spec}$")
plt.ylim(1e-6, 3)
plt.semilogy()

fig = plt.gcf()
fig.set_size_inches(5.52, 3.41)
fig.savefig("plots/out/ring_3.pdf", bbox_inches="tight")

plt.clf()

# fourth plot
max_sigma_inset = 200

dimfile = dim_dir + "ring_random_maxlength.dat"
print("Reading file {}".format(dimfile))
rand_conn, sigma, dim = np.loadtxt(dimfile, unpack=True)
rand_conn_vals = np.unique(rand_conn)

fig, ax1 = plt.subplots()
linecycler = cycle(lines)

for rc in rand_conn_vals:
    mask = rc == rand_conn
    these_sigmas = sigma[mask]
    ax1.plot(these_sigmas, dim[mask], label="${:0.0f}$".format(rc), ls=next(linecycler))

ax1.legend(loc="lower right")
ax1.set_xlabel("$\\sigma$")
ax1.set_ylabel("$d_{\\rm spec}$")

plt.axhline(2, c="tab:green", ls="--")
ax1.set_xlim([0, np.max(sigma)])
ax1.set_ylim([0.0, 1.2])

ax2 = plt.axes([0, 0, 1, 1])
# Manually set the position and relative size of the inset axes within ax1
ip = InsetPosition(ax1, [0.08, 0.1, 0.45, 0.45])
ax2.set_axes_locator(ip)
# ax2.set_xlabel("$\\sigma$")
# ax2.set_ylabel("$d_{\\rm spec}$")

linecycler = cycle(lines)

for rc in rand_conn_vals:
    mask = (rc == rand_conn) & (sigma < max_sigma_inset)
    these_sigmas = sigma[mask]
    ax2.plot(these_sigmas, dim[mask], ls=next(linecycler))

fig = plt.gcf()
fig.set_size_inches(5.52, 3.41)
fig.savefig("plots/out/ring_4.pdf", bbox_inches="tight")

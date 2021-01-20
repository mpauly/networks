from itertools import cycle

import matplotlib.pyplot as plt
import numpy as np

dimfile = "data/dim_ring_random.dat"

reference_delta = 0.2
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
        sigma[mask], dim[mask], label="${:d}$".format(round(nc)), ls=next(linecycler)
    )

plt.legend()
plt.title("$\\delta = {}$, varying number of connections".format(reference_delta))
plt.xlabel("$\\sigma$")
plt.ylabel("$d_{\\rm spec}$")
plt.xlim(0, 200)
plt.savefig("plots/out/dim_ring_1.pdf")

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
plt.savefig("plots/out/dim_ring_2.pdf")

plt.clf()

# third plot

dimfile = "data/dim_ring.dat"
print("Reading file {}".format(dimfile))
diffusion_const, sigma, dim = np.loadtxt(dimfile, unpack=True)
diffusion_const_values = np.unique(diffusion_const)

chain_length = 100

for dc in diffusion_const_values:
    mask = dc == diffusion_const
    these_sigmas = sigma[mask]
    line = plt.plot(these_sigmas, dim[mask], label="$\\delta={:.2f}$".format(dc))

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
plt.savefig("plots/out/dim_ring_3.pdf")

plt.clf()

# fourth plot

dimfile = "data/dim_ring_random_maxlength.dat"
print("Reading file {}".format(dimfile))
rand_conn, sigma, dim = np.loadtxt(dimfile, unpack=True)
rand_conn_vals = np.unique(rand_conn)

for rc in rand_conn_vals:
    mask = rc == rand_conn
    these_sigmas = sigma[mask]
    plt.plot(these_sigmas, dim[mask], label="conns=${:0.0f}$".format(rc))

plt.legend()
plt.xlabel("$\\sigma$")
plt.ylabel("$d_{\\rm spec}$")
plt.savefig("plots/out/dim_ring_4.pdf")

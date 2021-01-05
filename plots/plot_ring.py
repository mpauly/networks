import matplotlib.pyplot as plt
import numpy as np

dimfile = "data/dim_ring.dat"

reference_delta = 0.2
reference_conn = 200

print("Reading file {}".format(dimfile))
# 1D chain plot
nr_conn, diffusion_const, sigma, dim = np.loadtxt(dimfile, unpack=True)
nr_conn_values = np.unique(nr_conn)
diffusion_const_values = np.unique(diffusion_const)

# first plot
for nc in nr_conn_values:
    mask = (nc == nr_conn) & np.equal(diffusion_const, reference_delta)
    plt.plot(sigma[mask], dim[mask], label="${:d}$".format(round(nc)))

plt.legend()
plt.title("$\\delta = {}$".format(reference_delta))
plt.xlabel("$\\sigma$")
plt.ylabel("$d_{\\rm spec}$")
plt.xlim(0, 200)
plt.savefig("plots/out/dim_ring_1.pdf")

plt.clf()

# second plot
for dc in diffusion_const_values:
    mask = (dc == diffusion_const) & np.equal(nr_conn, reference_conn)
    plt.plot(sigma[mask], dim[mask], label="${:.2f}$".format(dc))

plt.legend()
plt.title("${}$ connections".format(reference_conn))
plt.xlabel("$\\sigma$")
plt.ylabel("$d_{\\rm spec}$")
plt.xlim(0, 200)
plt.savefig("plots/out/dim_ring_2.pdf")
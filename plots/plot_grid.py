import matplotlib.pyplot as plt
import numpy as np

dimfile = "data/dim_grid.dat"

print("Reading file {}".format(dimfile))
# 1D chain plot
diffusion_const, sigma, dim = np.loadtxt(dimfile, unpack=True)
diffusion_const_values = np.unique(diffusion_const)

# second plot
for dc in diffusion_const_values:
    mask = dc == diffusion_const
    plt.plot(sigma[mask], dim[mask], label="$\\delta={:.2f}$".format(dc))

plt.legend()
plt.xlabel("$\\sigma$")
plt.ylabel("$d_{\\rm spec}$")
# plt.xlim(0, 1000)
plt.savefig("plots/out/dim_grid.pdf")
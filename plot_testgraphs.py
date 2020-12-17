import matplotlib.pyplot as plt
import numpy as np

# 1D chain plot
sigma, dim = np.loadtxt('data/dim_1d_ring_26.dat', unpack=True)
mask = ~np.isnan(dim)
plt.plot(sigma[mask], dim[mask])
plt.axhline(1, ls='--')
plt.axvline(26, ls='--')
plt.axvline(52, ls='--')

plt.xlabel('$\\sigma$')
plt.ylabel('$d_{\\rm spec}$')
plt.savefig('plots/dim_1d_ring_26.png')

# 1D chain plot
plt.clf()
sigma, dim = np.loadtxt('data/dim_1d_ring_100.dat', unpack=True)
mask = ~np.isnan(dim)
plt.plot(sigma[mask], dim[mask])
sigma, dim = np.loadtxt('data/dim_1d_ring_100_with_random_50.dat', unpack=True)
mask = ~np.isnan(dim)
plt.plot(sigma[mask], dim[mask])
plt.axhline(1, ls='--')
plt.axvline(100, ls='--')

plt.xlabel('$\\sigma$')
plt.ylabel('$d_{\\rm spec}$')
plt.ylim([-0.2, 1.2])
plt.savefig('plots/dim_1d_ring_100_with_random_50.png')

# 2D plot
plt.clf()
sigma, dim = np.loadtxt('data/dim_2d_lattice_100.dat', unpack=True)

mask = ~np.isnan(dim)

plt.plot(sigma[mask], dim[mask])
plt.axhline(2, ls='--')

plt.xlabel('$\\sigma$')
plt.ylabel('$d_{\\rm spec}$')
plt.savefig('plots/dim_2d_lattice_100.png')

# 3D plot
plt.clf()
sigma, dim = np.loadtxt('data/dim_3d_lattice_100.dat', unpack=True)

mask = ~np.isnan(dim)

plt.plot(sigma[mask], dim[mask])
plt.axhline(3, ls='--')

plt.xlabel('$\\sigma$')
plt.ylabel('$d_{\\rm spec}$')
plt.savefig('plots/dim_3d_lattice_100.png')

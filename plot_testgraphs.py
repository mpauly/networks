import matplotlib.pyplot as plt
import numpy as np

# 1D chain plot
sigma, dim = np.loadtxt('data/dimension_1d_short.dat', unpack=True)
mask = ~np.isnan(dim)
plt.plot(sigma[mask], dim[mask])
plt.axhline(1, ls='--')
plt.axvline(26, ls='--')
plt.axvline(52, ls='--')

plt.xlabel('$\\sigma$')
plt.ylabel('$d_{\\rm spec}$')
plt.savefig('plots/dimension_1d_short.png')

# 1D chain plot
plt.clf()
sigma, dim = np.loadtxt('data/dimension_1d.dat', unpack=True)
mask = ~np.isnan(dim)
plt.plot(sigma[mask], dim[mask])
sigma, dim = np.loadtxt('data/dimension_1d_connected.dat', unpack=True)
mask = ~np.isnan(dim)
plt.plot(sigma[mask], dim[mask])
plt.axhline(1, ls='--')
plt.axvline(100, ls='--')

plt.xlabel('$\\sigma$')
plt.ylabel('$d_{\\rm spec}$')
plt.savefig('plots/dimension_1d.png')

# 2D plot
plt.clf()
sigma, dim = np.loadtxt('data/dimension_2d.dat', unpack=True)

mask = ~np.isnan(dim)

plt.plot(sigma[mask], dim[mask])
plt.axhline(2, ls='--')

plt.xlabel('$\\sigma$')
plt.ylabel('$d_{\\rm spec}$')
plt.savefig('plots/dimension_2d.png')

# 3D plot
plt.clf()
sigma, dim = np.loadtxt('data/dimension_3d.dat', unpack=True)

mask = ~np.isnan(dim)

plt.plot(sigma[mask], dim[mask])
plt.axhline(3, ls='--')

plt.xlabel('$\\sigma$')
plt.ylabel('$d_{\\rm spec}$')
plt.savefig('plots/dimension_3d.png')

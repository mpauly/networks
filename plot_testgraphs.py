import matplotlib.pyplot as plt
import numpy as np

sigma, dim = np.loadtxt('data/dimension_1d.dat', unpack=True)

mask = ~np.isnan(dim)

plt.plot(sigma[mask], dim[mask])
plt.axhline(1, ls='--')
plt.axvline(26, ls='--')
plt.axvline(52, ls='--')

plt.xlabel('$\\sigma$')
plt.ylabel('$d_{\\rm spec}$')
plt.savefig('plots/dimension_1d.png')

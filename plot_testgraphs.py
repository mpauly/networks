import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt('dimension_3D.dat')

line = data[0, 1::]
nan_mask = np.isfinite(line)
plt.plot(line[nan_mask])

plt.axhline(3, ls='--')

plt.xlabel('$\\sigma$')
plt.ylabel('$d_{\\rm spec}$')
plt.savefig('dimension_3d.png')

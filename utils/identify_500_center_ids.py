import numpy as np

data = np.loadtxt('data/other/nodes_cs_minkowski.csv')

centerdata = data[(data[:,1]>0.48) & (data[:,1]<0.52)&(data[:,2]>0.48)&(data[:,2]<0.52)]
ids = np.random.choice(centerdata[:,0], 500).astype(int)
np.savetxt('center_ids.txt', ids, fmt="%d")
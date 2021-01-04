import networkx as nx
import pickle

f = open("graphs/data/sub-0025864_ses-1_dwi_DS72784.gpickle", mode="rb")
outfile = open("graphs/data/brain_edgeslist.txt", mode="w")
g = pickle.load(f, encoding="bytes")
d = g.__dict__

nr_nodes = max(list(d[b"node"].keys()))

for i in range(nr_nodes):
    for j in range(i, nr_nodes):
        if j + 1 in d[b"edge"][i + 1]:
            outfile.write(
                "{} {} {}\n".format(i, j, d[b"edge"][i + 1][j + 1][b"weight"])
            )

f.close()
outfile.close()

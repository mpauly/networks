import glob
import os.path
import sys

import numpy as np

if len(sys.argv) != 2:
    raise ValueError("Please provide one directory full of walk files")

walk_dir = sys.argv[1]

outfile_name = walk_dir.replace("walks/", "data/dim_")
if outfile_name[-1] == "/":
    outfile_name = outfile_name[:-1]
outfile_name += ".dat"

outfile = open(outfile_name, "w")
outfile.write("# Merged dimension file - walks from {}\n".format(walk_dir))

walk_dir = os.path.join(walk_dir, "*.dat")

nr_of_walks = 0

for walk_file in glob.glob(walk_dir):
    with open(walk_file, "r") as f:
        l = f.readline()
        while l:
            if "start_node" in l:
                start_node = int(l.split("\t")[1])
                break
            l = f.readline()
        while l:
            if "dimension" in l:
                l = f.readline()
                dimension = np.fromstring(l, sep="\t")
                # trim the first and the last data point
                dimension = dimension[1:-2]
                break
            l = f.readline()
        print("Found walk {} with startnode {}".format(walk_file, start_node))
        # reformat numpy arrays
        for sigma, d in enumerate(np.nditer(dimension)):
            outfile.write("{}\t{}\t{}\n".format(start_node, sigma + 1, d))
        nr_of_walks += 1
        f.close()

outfile.close()

print("Merged a total of {} walks into file {}".format(nr_of_walks, outfile_name))

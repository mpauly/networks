infile = open(
    "data/raw/exported-traced-adjacencies-v1.2/traced-total-connections.csv", "r"
)
outfile = open(
    "data/raw/exported-traced-adjacencies-v1.2/traced-total-connections-renumbered.csv",
    "w",
)

infile.readline()

mapping = {}
counter = 0

for line in infile.readlines():
    inId, outId, weight = line.split(",")
    inId = mapping.get(inId, counter)
    if inId == counter:
        counter += 1
    outId = mapping.get(outId, counter)
    if outId == counter:
        counter += 1
    outfile.write("{},{},{}".format(inId, outId, weight))

infile.close()
outfile.close()

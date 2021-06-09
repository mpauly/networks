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
    if inId not in mapping:
        mapping[inId] = counter
        counter += 1
    inId = mapping[inId]
    if outId not in mapping:
        mapping[outId] = counter
        counter += 1
    outId = mapping[outId]
    outfile.write("{},{},{}".format(inId, outId, weight))

infile.close()
outfile.close()

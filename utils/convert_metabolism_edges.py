from itertools import product

import pandas as pd
from lxml import objectify

xml_data = objectify.parse("graphs/data/metabolism.kgml")

cmpds = xml_data.xpath("/pathway/entry[@type='compound']")
reactions = xml_data.xpath("/pathway/reaction")

edgecount = 0

nodes = set()

SEPARATOR = ","

outfile = open("graphs/data/metabolism.csv", "w")

for r in reactions:
    substrates = [s.get("id") for s in r.substrate]
    products = [p.get("id") for p in r.product]
    for edge in product(substrates, products):
        outfile.write(edge[0])
        outfile.write(SEPARATOR)
        outfile.write(edge[1])
        outfile.write("\n")
        edgecount += 1
        nodes.add(edge[0])
        nodes.add(edge[1])

print("total of {} nodes and {} edges".format(len(nodes), edgecount))

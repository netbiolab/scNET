#!/usr/bin/python3 -u
# usage: ./code.py [DIRECTORY OF RESULT FILE(first)] [diretory/of/BSorSV_fit/*_fit]
#./integration_1.py joined_net ./*_fit
#sorted to have a bigger lls value on the first column

import sys
print('modified code!!')
pairdata = dict()
for nw_index in range(len(sys.argv[2:])):
	with open(sys.argv[2:][nw_index]) as NW:
		for line in NW:
			a, b, lls = line.strip().split("\t")
			pair = tuple(sorted([a,b]))
			if pair not in pairdata:
				pairdata[pair] = [None] * len(sys.argv[2:])
			pairdata[pair][nw_index] = lls

output = open(sys.argv[1], "w")
for k in pairdata:
	print(k[0], k[1], "\t".join([x if x else "NA" for x in pairdata[k]]), sep="\t", file=output)

output.close()

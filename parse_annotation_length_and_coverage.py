#!/usr/bin/env python

import sys

# filtered gff3 from stdin
# counts of coverage vs length to stdout

count_per_length = {}
for l in sys.stdin:
	parts = l.strip().split("\t")
	name = parts[0]
	coverage = int(name.split("_")[1].replace("coverage", ""))
	length = int(parts[4]) - int(parts[3])
	if length not in count_per_length: count_per_length[length] = 0
	count_per_length[length] += coverage

lengths = [l for l in count_per_length.keys()]
lengths.sort()

#print("length,coverage")
for length in lengths:
	for i in range(0, count_per_length[length]):
		print(length)
	#print(str(length) + "," + str(count_per_length[length]))

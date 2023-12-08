#!/usr/bin/env python

import sys
import random

ground_truth_counts_file = sys.argv[1]
# fasta from stdin
# fasta to stdout

def get_random_seq(l):
	result = "".join("ATCG"[random.randint(0, 3)] for i in range(0, l))
	return result

ground_truth_counts = []
result_seq = ""
current_name = ""
for l in sys.stdin:
	if l[0] == ">":
		current_name = l[1:].strip()
	else:
		seq = l.strip()
		copy_count = random.randint(1, 30)
		for i in range(0, copy_count):
			result_seq += seq
		ground_truth_counts.append((current_name, copy_count))

with open(ground_truth_counts_file, "w") as f:
	for count in ground_truth_counts:
		f.write(count[0] + "\t" + str(count[1]) + "\n")

print(">concatenated_sequence")
print(get_random_seq(100000) + result_seq + get_random_seq(100000))

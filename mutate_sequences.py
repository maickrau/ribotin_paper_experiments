#!/usr/bin/env python

import sys
import random

substitution_rate = float(sys.argv[1]) * 4.0/3.0 # multiply because it picks a random base, so 1/4 probability of doing nothing
insertion_rate = float(sys.argv[2]) / 5.0 # divide by 5 because length is random from 1 to 10
deletion_rate = float(sys.argv[3]) / 5.0 # divide by 5 because length is random from 1 to 10
# fasta from stdin
# fasta to stdout

def mutate(raw_seq):
	result = ""
	pos = 0
	while pos < len(raw_seq):
		if random.uniform(0, 1) < substitution_rate:
			result += "ATCG"[random.randint(0, 3)]
			pos += 1
			continue
		if random.uniform(0, 1) < insertion_rate:
			insertion_length = random.randint(1,10)
			for i in range(0, insertion_length):
				result += "ATCG"[random.randint(0, 3)]
			continue
		if random.uniform(0, 1) < deletion_rate:
			pos += random.randint(0, 10)
			continue
		result += raw_seq[pos]
		pos += 1
	return result

sequence_number = 0
for l in sys.stdin:
	if l[0] == ">": continue
	seq = mutate(l.strip())
	print(">sequence" + str(sequence_number))
	print(seq)
	sequence_number += 1

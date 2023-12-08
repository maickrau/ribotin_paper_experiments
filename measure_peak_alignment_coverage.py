#!/usr/bin/env python

import sys

# paf from stdin
# result to stdout

min_alignment_length = 3000
reference_length = 7195

coverage = []
for i in range(0, reference_length):
	coverage.append(0)

positions = []
for l in sys.stdin:
	parts = l.strip().split("\t")
	read_start_pos = int(parts[2])
	read_end_pos = int(parts[3])
	ref_start_pos = int(parts[7])
	ref_end_pos = int(parts[8])
	if read_end_pos - read_start_pos < min_alignment_length: continue
	if ref_end_pos - ref_start_pos < min_alignment_length: continue
	for i in range(ref_start_pos, ref_end_pos):
		coverage[i] += 1

max_coverage = 0
for cov in coverage:
	if cov > max_coverage: max_coverage = cov

print(max_coverage)

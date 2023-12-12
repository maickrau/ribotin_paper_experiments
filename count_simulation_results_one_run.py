#!/usr/bin/env python

import sys

aln_file = sys.argv[1]
ground_truth_counts_file = sys.argv[2]
ribotin_morphs_file = sys.argv[3]
# summary to stdout

matches = set()
truth_has_match = set()
ribotin_has_match = set()
truth_has_multiple_match = set()
ribotin_has_multiple_match = set()
with open(aln_file) as f:
	for l in f:
		parts = l.strip().split("\t")
		ribotinname = parts[0]
		truthname = parts[5]
		ribotinlength = int(parts[1])
		ribotinstart = int(parts[2])
		ribotinend = int(parts[3])
		if ribotinend-ribotinstart < float(ribotinlength) * 0.99: continue
		truthlength = int(parts[6])
		truthstart = int(parts[7])
		truthend = int(parts[8])
		if truthend-truthstart < float(truthlength) * 0.99: continue
		mismatches = int(parts[12][5:])
		if ribotinstart > 0: mismatches += ribotinstart
		if ribotinend < ribotinlength: mismatches += ribotinlength-ribotinend
		if truthstart > 0: mismatches += truthstart
		if truthend < truthlength: mismatches += truthlength-truthend
		if mismatches > ribotinlength * 0.01: continue
		if mismatches > truthlength * 0.01: continue
		matches.add((truthname, ribotinname))
		if truthname in truth_has_match: truth_has_multiple_match.add(truthname)
		if ribotinname in ribotin_has_match: ribotin_has_multiple_match.add(ribotinname)
		truth_has_match.add(truthname)
		ribotin_has_match.add(ribotinname)

true_positives = len(truth_has_match)

print("new run")

count_truth_copycountfive = 0
true_positives_copycountfive = 0
false_negatives_copycountfive = 0
false_negatives = 0
truth_copy_count = {}
with open(ground_truth_counts_file) as f:
	for l in f:
		parts = l.strip().split("\t")
		truth_copy_count[parts[0]] = int(parts[1])
		if parts[0] not in truth_has_match:
			false_negatives += 1
			print("FN copy count " + str(int(parts[1])) + " " + parts[0])
		if int(parts[1]) >= 5:
			count_truth_copycountfive += 1
			if parts[0] in truth_has_match:
				true_positives_copycountfive += 1
			else:
				false_negatives_copycountfive += 1

false_positives_coveragethirty = 0
false_positives = 0
count_ribotin_coveragethirty = 0
ribotin_morph_coverage = {}
with open(ribotin_morphs_file) as f:
	for l in f:
		if l[0] == '>':
			name = l[1:].strip()
			coverage = int(name.split("_")[1].replace("coverage", ""))
			if name not in ribotin_has_match:
				false_positives += 1
				print("FP coverage " + str(coverage) + " " + name)
			ribotin_morph_coverage[name] = coverage
			if coverage >= 30:
				count_ribotin_coveragethirty += 1
				if name not in ribotin_has_match: false_positives_coveragethirty += 1

print("count ribotin morphs: " + str(len(ribotin_morph_coverage)))
print("sum TP: " + str(true_positives))
print("sum FN: " + str(false_negatives))
print("sum FP: " + str(false_positives))
print("count copycount>5: " + str(count_truth_copycountfive))
print("sum TP>5: " + str(true_positives_copycountfive))
print("sum FN>5: " + str(false_negatives_copycountfive))
print("count coverage>30: " + str(count_ribotin_coveragethirty))
print("sum FP>30: " + str(false_positives_coveragethirty))

with open(aln_file) as f:
	for l in f:
		parts = l.strip().split("\t")
		ribotinname = parts[0]
		truthname = parts[5]
		ribotinlength = int(parts[1])
		ribotinstart = int(parts[2])
		ribotinend = int(parts[3])
		if ribotinend-ribotinstart < float(ribotinlength) * 0.99: continue
		truthlength = int(parts[6])
		truthstart = int(parts[7])
		truthend = int(parts[8])
		if truthend-truthstart < float(truthlength) * 0.99: continue
		mismatches = int(parts[12][5:])
		if ribotinstart > 0: mismatches += ribotinstart
		if ribotinend < ribotinlength: mismatches += ribotinlength-ribotinend
		if truthstart > 0: mismatches += truthstart
		if truthend < truthlength: mismatches += truthlength-truthend
		if mismatches > ribotinlength * 0.01: continue
		if mismatches > truthlength * 0.01: continue
		if truthname in truth_has_multiple_match or ribotinname in ribotin_has_multiple_match:
			print("skip match " + truthname + " " + ribotinname + ", not unique")
			continue
		print("unique aln with num mismatches " + str(mismatches))
		print("unique aln with copy count " + str(truth_copy_count[truthname]) + " coverage " + str(ribotin_morph_coverage[ribotinname]))

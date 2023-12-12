#!/usr/bin/env python

import sys
import random

gfa_file = sys.argv[1]
# paths to stdout

def revcomp(s):
	comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
	return "".join(comp[c] for c in s[::-1])

def revnode(n):
	assert n[0] == ">" or n[1] == "<"
	return ("<" if n[0] == ">" else ">") + n[1:]

def get_seq(node_seqs, pos):
	if pos[0] == ">":
		return node_seqs[pos[1:]]
	else:
		assert pos[0] == "<"
		return revcomp(node_seqs[pos[1:]])

def sample_path(node_seqs, core_node_nodes, alleles):
	pathseq = ""
	pathstr = ""
	for j in range(0, len(core_node_nodes)):
		allele_index = random.randint(0, len(alleles[j])-1)
		for node in alleles[j][allele_index]:
			pathseq += get_seq(node_seqs, node)
			pathstr += node
		pathseq += get_seq(node_seqs, core_node_nodes[j])
		pathstr += core_node_nodes[j]
		assert len(alleles[j]) >= 1
	allele_index = random.randint(0, len(alleles[-1])-1)
	for node in alleles[-1][allele_index]:
		pathseq += get_seq(node_seqs, node)
		pathstr += node
	return pathseq

def get_string_edit_distance(seq1, seq2):
	if len(seq1) < len(seq2): return get_string_edit_distance(seq2, seq1)
	assert len(seq2) <= len(seq1)
	row = []
	for i in range(0, len(seq2)+1):
		row.append(i)
	for i in range(0, len(seq1)):
		next_row = []
		next_row.append(row[0]+1)
		for j in range(0, len(seq2)):
			min_score = min(min(row[j+1]+1, row[j]+(0 if seq1[i] == seq2[j] else 1)), next_row[-1]+1)
			next_row.append(min_score)
		row = next_row
	return row[-1]

def get_allele_edit_distance(node_seqs, allele1, allele2):
	seq1 = ""
	seq2 = ""
	for node in allele1:
		seq1 += get_seq(node_seqs, node)
	for node in allele2:
		seq2 += get_seq(node_seqs, node)
	return get_string_edit_distance(seq1, seq2)

node_seqs = {}
edges = {}
paths = []
with open(gfa_file) as f:
	for l in f:
		parts = l.strip().split("\t")
		if parts[0] == "S":
			node_seqs[parts[1]] = parts[2]
		elif parts[0] == "P":
			paths.append([(">" if c[-1] == "+" else "<") + c[:-1] for c in parts[2].split(",")])

not_unique_nodes = set()
not_core_nodes = set()
for path in paths:
	found_here = set()
	for node in path:
		if node[1:] in found_here: not_unique_nodes.add(node[1:])
		found_here.add(node[1:])
	for node in node_seqs:
		if node not in found_here: not_core_nodes.add(node)

core_nodes = set()
core_node_order = {}
core_node_nodes = []
for node in paths[0]:
	if node[1:] in not_core_nodes or node[1:] in not_unique_nodes: continue
	core_nodes.add(node[1:])
	core_node_order[node[1:]] = len(core_node_nodes)
	core_node_nodes.append(node)

alleles = []
for i in range(0, len(core_nodes)):
	alleles.append([])
alleles.append([])

for path in paths:
	last_core = -1
	for i in range(0, len(path)):
		if path[i][1:] not in core_nodes: continue
		assert last_core == -1 or core_node_order[path[last_core][1:]]+1 == core_node_order[path[i][1:]]
		alleles[core_node_order[path[i][1:]]].append(path[last_core+1:i])
		last_core = i
	alleles[-1].append(path[last_core+1:])

total_average_distance_sum = 0.0
for i in range(0, len(alleles)):
	distance_sum = 0
	for j in range(0, len(alleles[i])):
		for k in range(j+1, len(alleles[i])):
			distance_sum += get_allele_edit_distance(node_seqs, alleles[i][j], alleles[i][k]) * 2 # twice since it could be sampled either way
	average_distance_from_here = float(distance_sum)/float(len(alleles[i])*len(alleles[i]))
	total_average_distance_sum += average_distance_from_here
	print("site " + str(i) + " average edit distance " + str(average_distance_from_here))
print(total_average_distance_sum)

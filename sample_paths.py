#!/usr/bin/env python

import sys
import random

gfa_file = sys.argv[1]
num_paths = int(sys.argv[2])
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

for i in range(0, num_paths):
	pathseq = sample_path(node_seqs, core_node_nodes, alleles)
	print(">path" + str(i))
	print(pathseq)

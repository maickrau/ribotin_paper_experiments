#!/usr/bin/env python

import sys

gfa_file = sys.argv[1]
# count to stdout

def revnode(n):
	return (">" if n[0] == "<" else "<") + n[1:]

def getone(s):
	for c in s: return c

def find_biallelic_bubble(s, edges):
	if len(edges[s]) != 2: return None
	endnode = None
	for edge in edges[s]:
		if len(edges[revnode(edge)]) != 1: return None
		if len(edges[edge]) != 1: return None
		if endnode is None: endnode = getone(edges[edge])
		if endnode != getone(edges[edge]): return None
	if len(edges[revnode(endnode)]) != 2: return None
	return endnode

edges = {}
start_nodes = set()

with open(gfa_file) as f:
	for l in f:
		parts = l.strip().split("\t")
		if parts[0] == "S":
			start_nodes.add(">" + parts[1])
			start_nodes.add("<" + parts[1])
			if (">" + parts[1]) not in edges: edges[">" + parts[1]] = set()
			if ("<" + parts[1]) not in edges: edges["<" + parts[1]] = set()
		if parts[0] == "L":
			fromnode = (">" if parts[2] == "+" else "<") + parts[1]
			tonode = (">" if parts[4] == "+" else "<") + parts[3]
			if fromnode not in edges: edges[fromnode] = set()
			edges[fromnode].add(tonode)
			if revnode(tonode) not in edges: edges[revnode(tonode)] = set()
			edges[revnode(tonode)].add(revnode(fromnode))

count = 0
for node in start_nodes:
	bubble = find_biallelic_bubble(node, edges)
	if not bubble: continue
	count += 1

print("raw count: " + str(count))
print("bubble count: " + str(count/2))
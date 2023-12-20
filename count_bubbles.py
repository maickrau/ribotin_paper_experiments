#!/usr/bin/env python

import sys

gfa_file = sys.argv[1]
# count to stdout

def revnode(n):
	return (">" if n[0] == "<" else "<") + n[1:]

def getone(s):
	for c in s: return c

# Detecting Superbubbles in Assembly Graphs, Onodera et al 2013
# fig. 5
def find_bubble(s, edges):
	if s not in edges: return None
	if len(edges[s]) < 2: return None
	S = [s]
	visited = set()
	seen = set()
	seen.add(s)
	while len(S) > 0:
		v = S.pop()
		assert v in seen
		seen.remove(v)
		assert v not in visited
		visited.add(v)
		if v not in edges: return None
		if len(edges[v]) == 0: return None
		for u in edges[v]:
			if u[1:] == v[1:]: return None
			if revnode(u) in visited: return None
			if u == s: return None
			assert u not in visited
			seen.add(u)
			assert revnode(u) in edges
			assert len(edges[revnode(u)]) >= 1
			has_nonvisited_parent = False
			for parent_edge in edges[revnode(u)]:
				parent = revnode(parent_edge)
				if parent not in visited: has_nonvisited_parent = True
			if not has_nonvisited_parent: S.append(u)
		if len(S) == 1 and len(seen) == 1 and S[0] == getone(seen):
			t = S.pop()
			if t in edges:
				for edge in edges[t]:
					if edge == s: return None
			return (s, t)
	return None

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
	bubble = find_bubble(node, edges)
	if not bubble: continue
	count += 1

print("raw count: " + str(count))
print("bubble count: " + str(count/2))
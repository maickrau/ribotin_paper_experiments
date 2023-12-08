#!/usr/bin/env python

import sys

morph_fasta = sys.argv[1]
morph_annotation = sys.argv[2]
# fasta to stdout

read_sequences = {}
current_name = ""
current_seq = ""
with open(morph_fasta) as f:
	for l in f:
		if l[0] == ">":
			if len(current_seq) > 0: read_sequences[current_name] = current_seq
			current_name = l[1:].strip()
			current_seq = ""
		else:
			current_seq += l.strip()
if len(current_seq) > 0: read_sequences[current_name] = current_seq

with open(morph_annotation) as f:
	for l in f:
		parts = l.strip().split("\t")
		if len(parts) < 8: continue
		morphname = parts[0]
		startpos = int(parts[3])
		endpos = int(parts[4])
		tags = parts[8].split(";")
		annotation_id = ""
		annotation_product = ""
		annotation_note = ""
		for tag in tags:
			parts = tag.split("=")
			if parts[0] == "ID": annotation_id = tag.replace(' ', '_')
			if parts[0] == "Note": annotation_note = tag.replace(' ', '_')
			if parts[0] == "product": annotation_product = tag.replace(' ', '_')
		print(">" + morphname + "_" + str(startpos) + "_" + str(endpos) + "_" + annotation_id + "_" + annotation_note + "_" + annotation_product)
		print(read_sequences[morphname][startpos:endpos])

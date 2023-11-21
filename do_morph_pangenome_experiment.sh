#!/usr/bin/sh

mkdir hg002_morph_pangenome
cd hg002_morph_pangenome

cp ../hg002_experiment/out_ref/morphs.fa .
cp ../hg002_experiment/out_ref/consensus.fa .
cp ../hg002_experiment/out_ref/morph-annotations.gff3 .

# get the sequences of the 18s, 5.8s and 28s genes
grep -v '>' < ribotin_folder/template_seqs/rDNA_one_unit.fasta | tr -d '\n' | head -c 5526 | tail -c 1868 > seq_18s.fa
grep -v '>' < ribotin_folder/template_seqs/rDNA_one_unit.fasta | tr -d '\n' | head -c 6753 | tail -c 156 > seq_5.8s.fa
grep -v '>' < ribotin_folder/template_seqs/rDNA_one_unit.fasta | tr -d '\n' | head -c 12971 | tail -c 5050 > seq_28s.fa

minimap2 -t 8 -X -c -x asm5 morphs.fa morphs.fa > alns.paf
seqwish -s morphs.fa -p alns.paf -g graph.gfa
GraphAligner -x vg -g graph.gfa -f consensus.fa -a consensus_aln.gaf
GraphAligner -x vg -g graph.gfa -f seq_18s.fa -f seq_5.8s.fa -f seq_28s.fa -a gene_aln.gaf

awk '$3=="rRNA"' < morph-annotations.gff3 | grep "18S" | awk '{print $1 "\t" $5-$4 "bp";}' > gene_lengths_18s.txt
awk '$3=="rRNA"' < morph-annotations.gff3 | grep "5.8S" | awk '{print $1 "\t" $5-$4 "bp";}' > gene_lengths_5.8S.txt
awk '$3=="rRNA"' < morph-annotations.gff3 | grep "28S" | awk '{print $1 "\t" $5-$4 "bp";}' > gene_lengths_28S.txt

cd ..

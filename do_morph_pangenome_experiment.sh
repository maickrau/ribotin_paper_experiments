#!/usr/bin/sh

mkdir hg002_morph_pangenome
cd hg002_morph_pangenome

cp ../hg002_experiment/out_ref/morphs.fa .
cp ../hg002_experiment/out_ref/consensus.fa .
cp ../hg002_experiment/out_ref/morph-annotations.gff3 .

# get the sequences of the 18s, 5.8s and 28s genes
grep -v '>' < ../KY962518.1.fa | tr -d '\n' | head -c 5526 | tail -c 1868 | awk 'BEGIN{print ">KY962518.1_18s"}{print;}' > seq_18s.fa
grep -v '>' < ../KY962518.1.fa | tr -d '\n' | head -c 6753 | tail -c 156 | awk 'BEGIN{print ">KY962518.1_5.8s"}{print;}' > seq_5.8s.fa
grep -v '>' < ../KY962518.1.fa | tr -d '\n' | head -c 12971 | tail -c 5050 | awk 'BEGIN{print ">KY962518.1_28s"}{print;}' > seq_28s.fa

minimap2 -t 8 -X -c -x asm5 morphs.fa morphs.fa > alns.paf
seqwish -s morphs.fa -p alns.paf -g graph.gfa
GraphAligner -x vg -g graph.gfa -f consensus.fa -a consensus_aln.gaf
GraphAligner -x vg -g graph.gfa -f seq_18s.fa -f seq_5.8s.fa -f seq_28s.fa -a gene_aln.gaf
minimap2 -t 8 -c -x asm5 morphs.fa seq_18s.fa seq_5.8s.fa seq_28s.fa > alns_genes_to_morphs.paf
../extract_annotated_parts.py morphs.fa morph-annotations.gff3 > annotated_parts.fa
GraphAligner -x vg -g graph.gfa -f annotated_parts.fa -a annotation_alns.gaf
grep "exon" < morph-annotations.gff3 | grep "5.8S ribosomal RNA" > annotations_5.8s.gff3
grep "exon" < morph-annotations.gff3 | grep "18S ribosomal RNA" > annotations_18s.gff3
grep "exon" < morph-annotations.gff3 | grep "28S ribosomal RNA" > annotations_28s.gff3
../extract_annotated_parts.py morphs.fa annotations_5.8s.gff3 > morph_seqs_5.8s.fa
../extract_annotated_parts.py morphs.fa annotations_18s.gff3 > morph_seqs_18s.fa
../extract_annotated_parts.py morphs.fa annotations_28s.gff3 > morph_seqs_28s.fa

awk '$3=="rRNA"' < morph-annotations.gff3 | grep "18S" | awk '{print $1 "\t" $5-$4 "bp";}' > gene_lengths_18s.txt
awk '$3=="rRNA"' < morph-annotations.gff3 | grep "5.8S" | awk '{print $1 "\t" $5-$4 "bp";}' > gene_lengths_5.8S.txt
awk '$3=="rRNA"' < morph-annotations.gff3 | grep "28S" | awk '{print $1 "\t" $5-$4 "bp";}' > gene_lengths_28S.txt
grep -v '>' < morph_seqs_5.8s.fa | sort | uniq -c > gene_alleles_5.8s.txt
grep -v '>' < morph_seqs_18s.fa | sort | uniq -c > gene_alleles_18s.txt
grep -v '>' < morph_seqs_28s.fa | sort | uniq -c > gene_alleles_28s.txt

# number of nodes and edges: check with bandage

# variants per locus: count with bandage, get nodes with:
# grep "3'ETS" < annotation_alns.gaf | cut -f 6 | tr '<>' '\n' | sort | uniq | tr '\n' ',' | less
# and similarly for "5'ETS", "ITS-1", "ITS-2", "28S", "18S", "5.8S"

# longest stretch without variants
grep -P '^S' < graph.gfa | cut -f3 | wc -L

# number of bubbles
../count_bubbles.py graph.gfa
../count_biallelic_bubbles.py graph.gfa

# length histograms
awk '$3=="exon"' < morph-annotations.gff3 | grep "18S" | ../parse_annotation_length_and_coverage.py > length_histogram_18S.csv
awk '$3=="exon"' < morph-annotations.gff3 | grep "5.8S" | ../parse_annotation_length_and_coverage.py > length_histogram_5.8S.csv
awk '$3=="exon"' < morph-annotations.gff3 | grep "28S" | ../parse_annotation_length_and_coverage.py > length_histogram_28S.csv
awk '$3=="exon"' < morph-annotations.gff3 | grep "5'ETS" | ../parse_annotation_length_and_coverage.py > length_histogram_5primeETS.csv
awk '$3=="exon"' < morph-annotations.gff3 | grep "3'ETS" | ../parse_annotation_length_and_coverage.py > length_histogram_3primeETS.csv
awk '$3=="exon"' < morph-annotations.gff3 | grep "ITS-1" | ../parse_annotation_length_and_coverage.py > length_histogram_ITS-1.csv
awk '$3=="exon"' < morph-annotations.gff3 | grep "ITS-2" | ../parse_annotation_length_and_coverage.py > length_histogram_ITS-2.csv
awk '$3=="repeat_region"' < morph-annotations.gff3 | grep "SSR1" | ../parse_annotation_length_and_coverage.py > length_histogram_SSR1.csv
awk '$3=="repeat_region"' < morph-annotations.gff3 | grep "SSR2" | ../parse_annotation_length_and_coverage.py > length_histogram_SSR2.csv
awk '$3=="repeat_region"' < morph-annotations.gff3 | grep "SSR3" | ../parse_annotation_length_and_coverage.py > length_histogram_SSR3.csv
awk '$3=="tandem_repeat"' < morph-annotations.gff3 | grep "TR1" | ../parse_annotation_length_and_coverage.py > length_histogram_TR1.csv
awk '$3=="tandem_repeat"' < morph-annotations.gff3 | grep "TR2" | ../parse_annotation_length_and_coverage.py > length_histogram_TR2.csv
awk '$3=="repeat_region"' < morph-annotations.gff3 | grep "similar to Long Repeat 1" | ../parse_annotation_length_and_coverage.py > length_histogram_LR1.csv
awk '$3=="repeat_region"' < morph-annotations.gff3 | grep "LR2" | ../parse_annotation_length_and_coverage.py > length_histogram_LR2.csv

Rscript ../plot_morph_pangenome_length_histograms.Rscript

cd ..

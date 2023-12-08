#!/usr/bin/sh

mkdir hg002_morph_pangenome
cd hg002_morph_pangenome

cp ../hg002_experiment/out_ref/morphs.fa .
cp ../hg002_experiment/out_ref/consensus.fa .
cp ../hg002_experiment/out_ref/morph-annotations.gff3 .

# get the sequences of the 18s, 5.8s and 28s genes
grep -v '>' < ribotin_folder/template_seqs/rDNA_one_unit.fasta | tr -d '\n' | head -c 5526 | tail -c 1868 | awk 'BEGIN{print ">KY962518.1_18s"}{print;}' > seq_18s.fa
grep -v '>' < ribotin_folder/template_seqs/rDNA_one_unit.fasta | tr -d '\n' | head -c 6753 | tail -c 156 | awk 'BEGIN{print ">KY962518.1_5.8s"}{print;}' > seq_5.8s.fa
grep -v '>' < ribotin_folder/template_seqs/rDNA_one_unit.fasta | tr -d '\n' | head -c 12971 | tail -c 5050 | awk 'BEGIN{print ">KY962518.1_28s"}{print;}' > seq_28s.fa

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

cd ..

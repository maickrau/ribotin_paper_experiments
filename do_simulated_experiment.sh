#!/usr/bin/sh

mkdir simulated_experiment
cd simulated_experiment

minimap2 -t 4 -x asm5 -c --eqx ../chm13_major_morphs.fa ../chm13_major_morphs.fa > alns_ava.paf
seqwish -s ../chm13_major_morphs.fa -p alns_ava.paf -g chm13_morph_graph.gfa
../sample_paths.py chm13_morph_graph.fa 20 > paths.fa
../mutate_sequences.py 0.005 0.002 0.002 < paths.fa > mutated_paths.fa
../multiply_and_concatenate.py ground_truth_counts.txt < mutated_paths.fa > rDNA_sequence_ground_truth.fa
pbsim --strategy wgs --method qshmm --qshmm QSHMM-ONT-HQ.model --depth 30 --genome rDNA_sequence_ground_truth.fa --length-mean 80000 --length-sd 70000
mv sd_0001.fastq simulated_ont_reads.fq
pbsim --strategy wgs --method qshmm --qshmm QSHMM-RSII.model --depth 30 --genome rDNA_sequence_ground_truth.fa --pass-num 10 --length-mean 14000 --length-sd 2000
samtools view -b sd_0001.sam > sd_0001.bam
ccs sd_0001.bam simulated_hifi.bam

#!/usr/bin/sh

mkdir chm13_experiment
cd chm13_experiment

mkdir chm13_data
cd chm13_data
wget -O - https://sra-pub-src-2.s3.amazonaws.com/SRR11292120/m64062_190806_063919.fastq.1 | awk 'NR%4==1||NR%4==2' | tr '@' '>' | gzip >m64062_190806_063919.fa.gz &
wget -O - https://sra-pub-src-2.s3.amazonaws.com/SRR11292121/m64062_190803_042216.fastq.1 | awk 'NR%4==1||NR%4==2' | tr '@' '>' | gzip > m64062_190803_042216.fa.gz &
wget -O - https://sra-pub-src-2.s3.amazonaws.com/SRR11292122/m64062_190807_194840.fastq.1 | awk 'NR%4==1||NR%4==2' | tr '@' '>' | gzip > m64062_190807_194840.fa.gz &
wget -O - https://sra-pub-src-2.s3.amazonaws.com/SRR11292123/m64062_190804_172951.fastq.1 | awk 'NR%4==1||NR%4==2' | tr '@' '>' | gzip > m64062_190804_172951.fa.gz &
wget -O - https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/nanopore/rel8-guppy-5.0.7/reads.fastq.gz | zcat | awk 'NR%4==1||NR%4==2' | tr '@' '>' | /projappl/project_2006830/filter_fasta_by_length.py 5000 | gzip > ont_reads_5k.fa.gz &
wait
cd ..

mkdir verkko_asm
verkko -d verkko_asm --hifi chm13_data/m64062_190803_042216.fa.gz --hifi chm13_data/m64062_190804_172951.fa.gz --hifi chm13_data/m64062_190806_063919.fa.gz --hifi chm13_data/m64062_190807_194840.fa.gz --nano chm13_data/ont_reads_5k.fa.gz --slurm --snakeopts "--jobs 50"

ribotin-ref -x human -t 8 -i chm13_data/m64062_190803_042216.fa.gz -i chm13_data/m64062_190804_172951.fa.gz -i chm13_data/m64062_190806_063919.fa.gz -i chm13_data/m64062_190807_194840.fa.gz --nano chm13_data/ont_reads_5k.fa.gz -o out_ref
ribotin-verkko -x human -t 8 -i verkko_asm -o out_verkko_automatic

# manually pick the five rDNA clusters into cluster1.txt cluster2.txt cluster3.txt cluster4.txt cluster5.txt
ribotin-verkko -x human -t 8 -i verkko_asm -o out_verkko_manual -c cluster1.txt -c cluster2.txt -c cluster3.txt -c cluster4.txt -c cluster5.txt

mkdir alignments
minimap2 --eqx -x asm5 -c -t 8 chm13_major_morphs.fa out_ref/morphs.fa > alns_ref_t2t.paf
minimap2 --eqx -x asm5 -c -t 8 chm13_major_morphs.fa out_verkko_automatic*/morphs.fa > alns_verkko_t2t.paf
minimap2 --eqx -x asm5 -c -t 8 chm13_major_morphs.fa out_verkko_manual*/morphs.fa > alns_verkkomanual_t2t.paf
minimap2 --eqx -x asm5 -c -t 8 out_ref/morphs.fa out_verkko_automatic*/morphs.fa > alns_verkko_ref.paf
minimap2 --eqx -x asm5 -c -t 8 out_ref/consensus.fa out_ref/morphs.fa > alns_ref_consensus.paf
minimap2 --eqx -x asm5 -c -t 8 out_verkko_manual0/consensus.fa out_verkko_manual0/morphs.fa > alns_manual0_consensus0.paf
minimap2 --eqx -x asm5 -c -t 8 out_verkko_manual1/consensus.fa out_verkko_manual1/morphs.fa > alns_manual0_consensus1.paf
minimap2 --eqx -x asm5 -c -t 8 out_verkko_manual2/consensus.fa out_verkko_manual2/morphs.fa > alns_manual0_consensus2.paf
minimap2 --eqx -x asm5 -c -t 8 out_verkko_manual3/consensus.fa out_verkko_manual3/morphs.fa > alns_manual0_consensus3.paf
minimap2 --eqx -x asm5 -c -t 8 out_verkko_manual4/consensus.fa out_verkko_manual4/morphs.fa > alns_manual0_consensus4.paf

cd ..

#!/usr/bin/sh

mkdir arabidopsis_experiment
cd arabidopsis_experiment

mkdir arabidopsis_data
cd arabidopsis_data
wget -O - https://download.cncb.ac.cn/gsa/CRA004538/CRR302668/CRR302668.fastq.gz | gunzip | awk 'NR%4==1||NR%4==2' | tr '@' '>' | gzip > CRR302668.fa.gz &
wget -O - https://download.cncb.ac.cn/gsa/CRA004538/CRR302667/CRR302667.fastq.gz | gunzip | awk 'NR%4==1||NR%4==2' | tr '@' '>' | gzip > CRR302667.fa.gz &
wait


mkdir verkko_asm
verkko -d verkko_asm --hifi arabidopsis_data/CRR302668.fa.gz --nano arabidopsis_data/CRR302667.fa.gz --slurm --snakeopts "--jobs 50"

# first results using a human reference
ribotin-ref -x human -t 8 -i arabidopsis_data/CRR302668.fa.gz --nano arabidopsis_data/CRR302667.fa.gz -o out_ref_humanref
ribotin-verkko -x human -t 8 -i verkko_asm -o out_verkko_automatic_humanref

# then with a species-specific reference
# make the arabidopsis rDNA reference
MBG -t 8 -i arabidopsis_data/CRR302668.fa.gz -o mbg-graph.gfa -k 1001 -w 100 -a 1 -u 2 --error-masking=msat -r 15000 -R 4000
grep -P '^S' < mbg-graph.gfa | awk '{print ">" $2; print $3;}' > mbg-contigs.fa
# manually pick the arabidopsis rDNA contigs, put them in arabidopsis_rdna_kmers.fa
ribotin-ref -r arabidopsis_rdna_kmers.fa --approx-morphsize 10000 -t 8 -i arabidopsis_data/CRR302668.fa.gz --nano arabidopsis_data/CRR302667.fa.gz -o out_ref
ribotin-verkko --guess-tangles-using-reference arabidopsis_rdna_kmers.fa --approx-morphsize 10000 -t 8 -i verkko_asm -o out_verkko_automatic

# use the hifi reads as "ont" reads as well
ribotin-ref -r arabidopsis_rdna_kmers.fa --approx-morphsize 10000 -t 8 -i arabidopsis_data/CRR302668.fa.gz --nano arabidopsis_data/CRR302668.fa.gz -o out_ref_hifionly --morph-cluster-maxedit 10 --morph-recluster-minedit 1

mkdir alignments
minimap2 --eqx -x asm5 -c -t 8 out_ref/morphs.fa out_verkko_automatic*/morphs.fa > alns_verkko_ref.paf
minimap2 --eqx -x asm5 -c -t 8 out_ref/consensus.fa out_ref/morphs.fa > alns_ref_consensus.paf
minimap2 --eqx -x asm5 -c -t 8 out_ref_hifionly/consensus.fa out_ref_hifionly/morphs.fa > alns_ref_consensus_hifionly.paf

cd ..

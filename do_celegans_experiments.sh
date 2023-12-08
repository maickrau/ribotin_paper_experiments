#!/usr/bin/sh

mkdir celegans_experiment
cd celegans_experiment

mkdir celegans_data
cd celegans_data
wget -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR221/022/SRR22137522/SRR22137522_subreads.fastq.gz | gunzip | awk 'NR%4==1||NR%4==2' | tr '@' '>' | gzip > SRR22137522.fa.gz &
wget -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR221/023/SRR22137523/SRR22137523_subreads.fastq.gz | gunzip | awk 'NR%4==1||NR%4==2' | tr '@' '>' | gzip > SRR22137523.fa.gz &
wait


# first results using a human reference
ribotin-ref -x human -t 8 -i celegans_data/SRR22137522.fa.gz --nano celegans_data/SRR22137522.fa.gz -o out_ref_alt2_humanref
ribotin-ref -x human -t 8 -i celegans_data/SRR22137523.fa.gz --nano celegans_data/SRR22137523.fa.gz -o out_ref_alt1_humanref

# then with a species-specific reference
# make the c. elegans rDNA reference
MBG -t 8 -i celegans_data/SRR22137522.fa.gz -i celegans_data/SRR22137523.fa.gz -o mbg-graph.gfa -k 1001 -w 100 -a 1 -u 2 --error-masking=msat -r 15000 -R 4000
MBG -t 8 -i celegans_data/SRR22137522.fa.gz -o mbg-graph-alt2.gfa -k 1001 -w 100 -a 1 -u 2 --error-masking=msat -r 15000 -R 4000
MBG -t 8 -i celegans_data/SRR22137523.fa.gz -o mbg-graph-alt1.gfa -k 1001 -w 100 -a 1 -u 2 --error-masking=msat -r 15000 -R 4000
grep -P '^S' < mbg-graph.gfa | awk '{print ">" $2; print $3;}' > mbg-contigs.fa
# manually pick the c. elegans rDNA contigs, put them in celegans_rdna_kmers.fa
ribotin-ref -r celegans_rdna_kmers.fa --approx-morphsize 10000 --morph-cluster-maxedit 10 --morph-recluster-minedit 1 -t 8 -i celegans_data/SRR22137523.fa.gz --nano celegans_data/SRR22137523.fa.gz -o out_ref_alt1
ribotin-ref -r celegans_rdna_kmers.fa --approx-morphsize 10000 --morph-cluster-maxedit 10 --morph-recluster-minedit 1 -t 8 -i celegans_data/SRR22137522.fa.gz --nano celegans_data/SRR22137522.fa.gz -o out_ref_alt2

mkdir alignments
cd alignments
minimap2 --eqx -x asm5 -c -t 8 ../out_ref_alt1/consensus.fa out_ref_alt1/morphs.fa > alns_ref_consensus_alt1.paf
minimap2 --eqx -x asm5 -c -t 8 ../out_ref_alt2/consensus.fa out_ref_alt2/morphs.fa > alns_ref_consensus_alt2.paf
minimap2 --eqx -x asm5 -c -t 8 ../out_ref_alt1/consensus.fa out_ref_alt2/consensus.fa > alns_consensus_to_consensus.paf
minimap2 --eqx -x asm5 -c -t 8 ../out_ref_alt1/consensus.fa out_ref_alt2/morphs.fa > alns_ref_alt2_to_consensus_alt1.paf

#estimate copy counts by aligning all reads to consensus
minimap2 --eqx -x asm5 -c -t 24 ../out_ref_alt1/consensus.fa ../celegans_data/SRR22137523.fa.gz > alns_alt1reads_to_alt1consensus.paf
minimap2 --eqx -x asm5 -c -t 24 ../out_ref_alt1/consensus.fa ../celegans_data/SRR22137522.fa.gz > alns_alt2reads_to_alt1consensus.paf

cd ..

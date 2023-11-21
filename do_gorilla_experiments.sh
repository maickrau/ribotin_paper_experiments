#!/usr/bin/sh

mkdir gorilla_experiment
cd gorilla_experiment

mkdir gorilla_data
cd gorilla_data
wget -O - https://genomeark.s3.amazonaws.com/species/Gorilla_gorilla/mGorGor1/genomic_data/pacbio_hifi/m54329U_210319_174352.ccs.Q20.fastq.gz | gunzip | awk 'NR%4==1||NR%4==2' | tr '@' '>' | gzip > m54329U_210319_174352.fa.gz &
wget -O - https://genomeark.s3.amazonaws.com/species/Gorilla_gorilla/mGorGor1/genomic_data/pacbio_hifi/m64076_210208_175256.Q20.fastq.gz | gunzip | awk 'NR%4==1||NR%4==2' | tr '@' '>' | gzip > m64076_210208_175256.fa.gz &
wget -O - https://genomeark.s3.amazonaws.com/species/Gorilla_gorilla/mGorGor1/genomic_data/pacbio_hifi/m64076_210213_010909.Q20.fastq.gz | gunzip | awk 'NR%4==1||NR%4==2' | tr '@' '>' | gzip > m64076_210213_010909.fa.gz &
wget -O - https://genomeark.s3.amazonaws.com/species/Gorilla_gorilla/mGorGor1/genomic_data/ont/PAG67391_guppy-6.3.7-sup.fastq.gz | gunzip | awk 'NR%4==1||NR%4==2' | tr '@' '>' | gzip > PAG67391_guppy-6.3.7-sup.fa.gz &
wget -O - https://genomeark.s3.amazonaws.com/species/Gorilla_gorilla/mGorGor1/genomic_data/ont/PAG68076_guppy-6.3.7-sup.fastq.gz | gunzip | awk 'NR%4==1||NR%4==2' | tr '@' '>' | gzip > PAG68076_guppy-6.3.7-sup.fa.gz &
wait
cd ..

verkko -d verkko_asm --hifi gorilla_data/m54329U_210319_174352.fa.gz --hifi gorilla_data/m64076_210208_175256.fa.gz --hifi gorilla_data/m64076_210213_010909.fa.gz --nano gorilla_data/PAG67391_guppy-6.3.7-sup.fa.gz --nano gorilla_data/PAG68076_guppy-6.3.7-sup.fa.gz --slurm --snakeopts "--jobs 50"

# first results using a human reference
ribotin-ref -x human -t 8 -i gorilla_data/m54329U_210319_174352.fa.gz -i gorilla_data/m64076_210208_175256.fa.gz -i gorilla_data/m64076_210213_010909.fa.gz --nano gorilla_data/PAG67391_guppy-6.3.7-sup.fa.gz --nano gorilla_data/PAG68076_guppy-6.3.7-sup.fa.gz -o out_ref_humanref
ribotin-verkko -x human -t 8 -i verkko_asm -o out_verkko_automatic_humanref

# then with a species-specific reference
# make the gorilla rDNA reference
MBG -t 8 -i gorilla_data/m54329U_210319_174352.fa.gz -i gorilla_data/m64076_210208_175256.fa.gz -i gorilla_data/m64076_210213_010909.fa.gz -o mbg-graph.gfa -k 1001 -w 100 -a 1 -u 2 --error-masking=msat -r 15000 -R 4000
grep -P '^S' < mbg-graph.gfa | awk '{print ">" $2; print $3;}' > mbg-contigs.fa
# manually pick the gorilla rDNA contigs, put them in gorilla_rdna_kmers.fa
ribotin-ref -r gorilla_rdna_kmers.fa --approx-morphsize 45000 -t 8 -i gorilla_data/m54329U_210319_174352.fa.gz -i gorilla_data/m64076_210208_175256.fa.gz -i gorilla_data/m64076_210213_010909.fa.gz --nano gorilla_data/PAG67391_guppy-6.3.7-sup.fa.gz --nano gorilla_data/PAG68076_guppy-6.3.7-sup.fa.gz -o out_ref
ribotin-verkko -r gorilla_rdna_kmers.fa --approx-morphsize 45000 -t 8 -i verkko_asm -o out_verkko_automatic

mkdir alignments
minimap2 --eqx -x asm5 -c -t 8 out_ref/morphs.fa out_verkko_automatic*/morphs.fa > alns_verkko_ref.paf
minimap2 --eqx -x asm5 -c -t 8 out_ref/consensus.fa out_ref/morphs.fa > alns_ref_consensus.paf

cd ..

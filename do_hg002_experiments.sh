#!/usr/bin/sh

mkdir hg002_experiment
cd hg002_experiment

mkdir hg002_data
cd hg002_data

wget -O - https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/scratch/HG002/sequencing/hifi/m64011_190830_220126.Q20.fastq.gz | gunzip | awk 'NR%4==1||NR%4==2' | tr '@' '>' | gzip > m64011_190830_220126.Q20.fa.gz &
wget -O - https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/scratch/HG002/sequencing/hifi/m64011_190901_095311.Q20.fastq.gz | gunzip | awk 'NR%4==1||NR%4==2' | tr '@' '>' | gzip > m64011_190901_095311.Q20.fa.gz &
wget -O - https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/scratch/HG002/sequencing/hifi/m64012_190920_173625.Q20.fastq.gz | gunzip | awk 'NR%4==1||NR%4==2' | tr '@' '>' | gzip > m64012_190920_173625.Q20.fa.gz &
wget -O - https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/scratch/HG002/sequencing/ont/03_08_22_R941_HG002_1_Guppy_6.0.6_prom_sup.fastq.gz | gunzip | awk 'NR%4==1||NR%4==2' | tr '@' '>' | gzip > 03_08_22_R941_HG002_1_Guppy_6.0.6_prom_sup.fa.gz &
wget -O - https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/scratch/HG002/sequencing/ont/03_08_22_R941_HG002_2_Guppy_6.0.6_prom_sup.fastq.gz | gunzip | awk 'NR%4==1||NR%4==2' | tr '@' '>' | gzip > 03_08_22_R941_HG002_2_Guppy_6.0.6_prom_sup.fa.gz &
wait

# estimate ONT and hifi read length distributions
zcat 03_08_22_R941_HG002_2_Guppy_6.0.6_prom_sup.fa.gz | grep -v '>' | awk '{print length($0);}' > readlens_03_08_22_R941_HG002_2_Guppy_6.0.6_prom_sup.txt
awk '{sum += $1; count += 1;}END{print sum/count;}' < readlens_03_08_22_R941_HG002_2_Guppy_6.0.6_prom_sup.txt
# average ONT length 37431.9
awk '{sum += sqrt(($1-37431.9)*($1-37431.9)); count += 1;}END{print sum/count;}' < readlens_03_08_22_R941_HG002_2_Guppy_6.0.6_prom_sup.txt
# ONT standard deviation 32111.1
zcat m64011_190830_220126.Q20.fa.gz | grep -v '>' | awk '{print length($0);}' > readlens_m64011_190830_220126.Q20.txt
awk '{sum += $1; count += 1;}END{print sum/count;}' < readlens_m64011_190830_220126.Q20.txt
# average HiFi length 18545.4
awk '{sum += sqrt(($1-18545.4)*($1-18545.4)); count += 1;}END{print sum/count;}' < readlens_m64011_190830_220126.Q20.txt
# HiFi standard deviation 1622.27


cd ..

verkko -d verkko_asm --hifi hg002_data/m64011_190830_220126.Q20.fa.gz --hifi hg002_data/m64011_190901_095311.Q20.fa.gz --hifi hg002_data/m64012_190920_173625.Q20.fa.gz --nano hg002_data/03_08_22_R941_HG002_1_Guppy_6.0.6_prom_sup.fa.gz --nano hg002_data/03_08_22_R941_HG002_2_Guppy_6.0.6_prom_sup.fa.gz --slurm --snakeopts "--jobs 50"

ribotin-ref -x human -t 8 -i hg002_data/m64011_190830_220126.Q20.fa.gz -i hg002_data/m64011_190901_095311.Q20.fa.gz -i hg002_data/m64012_190920_173625.Q20.fa.gz --nano hg002_data/03_08_22_R941_HG002_1_Guppy_6.0.6_prom_sup.fa.gz --nano hg002_data/03_08_22_R941_HG002_2_Guppy_6.0.6_prom_sup.fa.gz -o out_ref
ribotin-verkko -x human -t 8 -i verkko_asm -o out_verkko_automatic

ribotin-ref -x human -t 8 -i hg002_data/m64011_190830_220126.Q20.fa.gz --nano hg002_data/03_08_22_R941_HG002_1_Guppy_6.0.6_prom_sup.fa.gz -o out_ref_lowcoverage

mkdir alignments
minimap2 --eqx -x asm5 -c -t 8 out_ref/morphs.fa out_verkko_automatic*/morphs.fa > alignments/alns_verkko_ref.paf
minimap2 --eqx -x asm5 -c -t 8 out_ref/consensus.fa out_ref/morphs.fa > alignments/alns_ref_consensus.paf
minimap2 --eqx -x asm5 -c -t 8 out_ref/morphs.fa out_ref_lowcoverage/morphs.fa > alignments/alns_lowcoverage_to_fullcoverage.paf

cd ..

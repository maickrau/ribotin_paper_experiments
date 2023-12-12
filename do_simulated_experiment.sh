#!/usr/bin/sh

mkdir simulated_experiment
cd simulated_experiment

cat ../chm13_major_morphs.fa ../KY962518.1.fa > rDNA_morphs.fa
minimap2 -t 4 -x asm5 -c --eqx rDNA_morphs.fa rDNA_morphs.fa > alns_ava.paf
seqwish -s rDNA_morphs.fa -p alns_ava.paf -g morph_graph.gfa
# check average divergence of mosaic sequences
../measure_sampled_path_average_empirical_edit_distance.py morph_graph.gfa
# average 1889.9 edits, assume 45000 average size -> 4.2% average divergence

mkdir run1
cd run1
../../sample_paths.py ../morph_graph.gfa 20 > paths.fa
../../mutate_sequences.py 0.001 0.0005 0.0005 < paths.fa > mutated_paths.fa
../../multiply_and_concatenate.py ground_truth_counts.txt < mutated_paths.fa > rDNA_sequence_ground_truth.fa
pbsim --strategy wgs --method qshmm --qshmm ../QSHMM-ONT-HQ.model --depth 30 --genome rDNA_sequence_ground_truth.fa --length-mean 40000 --length-sd 30000
mv sd_0001.fastq simulated_ont_reads.fq
# using the hifi lengths from do_hg002_experiment.sh results in a segmentation fault
# so round the lengths down to 15kbp
pbsim --strategy wgs --method qshmm --qshmm ../QSHMM-RSII.model --depth 10 --genome rDNA_sequence_ground_truth.fa --pass-num 10 --length-mean 15000 --length-sd 2000
samtools view -b sd_0001.sam > sd_0001.bam
ccs -j 8 sd_0001.bam simulated_hifi.bam
samtools fasta simulated_hifi.bam > simulated_hifi.fa

ribotin-ref -t 8 -x human -i simulated_hifi.fa --nano simulated_ont_reads.fq -o out
minimap2 --eqx -x asm5 -c -t 8 mutated_paths.fa out/morphs.fa > alns_ribotinref_to_ground_truth.paf
../../count_simulation_results_one_run.py alns_ribotinref_to_ground_truth.paf ground_truth_counts.txt out/morphs.fa > result_summary.txt
cat result_summary.txt >> ../all_runs_summary.txt

# repeat for multiple runs
cd ..

# average coverage per copy count
grep "unique aln" < all_runs_summary.txt | grep "copy count" | awk '{sum += $8; count += $6;}END{print sum/count;}'
# correlation of coverage vs copy count
grep "unique aln" < all_runs_summary.txt | grep "copy count" | cut -f 6,8 -d ' ' | awk '{count += 1; x += $1; y += $2; xx += $1*$1; yy += $2*$2; xy += $1*$2;}END{print ((xy/count)-(x/count*y/count))/(sqrt(xx/count-x*x/count/count)*sqrt(yy/count-y*y/count/count));}'
# average number of mismatches
grep "unique aln" < all_runs_summary.txt | grep "mismatches" | awk '{sum += $6; count += 1;}END{print sum/count;}'

# sensitivity and specificity stats
grep "new run" < all_runs_summary.txt | awk '{sum += 20;}END{print "total truth morphs: " sum;}' > final_summary_stats.txt
grep "count ribotin morphs" < all_runs_summary.txt | awk '{sum += $4;}END{print "total ribotin morphs: " sum;}' >> final_summary_stats.txt
grep "sum FP:" < all_runs_summary.txt | awk '{sum += $3;}END{print "sum FP: " sum;}' >> final_summary_stats.txt
grep "sum FN:" < all_runs_summary.txt | awk '{sum += $3;}END{print "sum FN: " sum;}' >> final_summary_stats.txt
grep "sum TP:" < all_runs_summary.txt | awk '{sum += $3;}END{print "sum TP: " sum;}' >> final_summary_stats.txt
grep "count copycount>5" < all_runs_summary.txt | awk '{sum += $3;}END{print "total copycount>5: " sum;}' >> final_summary_stats.txt
grep "sum FN>5:" < all_runs_summary.txt | awk '{sum += $3;}END{print "sum FN>5: " sum;}' >> final_summary_stats.txt
grep "sum TP>5:" < all_runs_summary.txt | awk '{sum += $3;}END{print "sum TP>5: " sum;}' >> final_summary_stats.txt
grep "count coverage>30" < all_runs_summary.txt | awk '{sum += $3;}END{print "total coverage>30: " sum;}' >> final_summary_stats.txt
grep "sum FP>30:" < all_runs_summary.txt | awk '{sum += $3;}END{print "sum FP>30: " sum;}' >> final_summary_stats.txt
awk '/total truth morphs/{truth_morphs=int($4)}/total ribotin morphs/{ribotin_morphs=int($4)}/sum FP:/{FP = $3}/sum FN:/{FN = $3}/sum TP:/{TP = $3}/total copycount>5/{truth_5 = int($3);}/total coverage>30/{ribotin_30 = $3}/sum FN>5/{FN5 = $3}/sum TP>5/{TP5 = $3}/sum FP>30/{FP30 = $3}END{print "sensitivity " TP/(TP+FN); print "specificity " TP/(TP+FP); print "sensitivity copycount5 " TP5/(TP5+FN5); print "specificity coverage30 " TP5/(TP5+FP30)}' < final_summary_stats.txt

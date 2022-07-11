#
for file_name in *.fastq;
do
fastqc -t 10 ${file_name} -o /dcl01/hansen/data/kabuki_mice/HiSeq325_ATAC_new/ATAC/QCreport;
done;

for file_name in *R1_001.fastq;
do
echo mapping ${file_name};
bowtie2 -p 15 -x /dcl01/hansen/data/bjornsson_ht22_chipseq_miseq/ks_chipseq/mm10/mm10 -1 ${file_name} -2 ${file_name/R1_001.fastq/}R2_001.fastq -S ${file_name/R1_001.fastq/}.sam;
done;



for file_name in *.sam;
do
samtools view -S -b ${file_name} > ${file_name/_.sam/}.bam;
done;

for file_name in *.bam;
do
echo sorting ${file_name};
samtools sort -@ 30 -m 2000M ${file_name} ${file_name/_.bam/}_sorted.bam;
done;

for file_name in *.bam;
do
echo sorting ${file_name};
samtools sort -@ 30 -m 2000M -o ${file_name/_.bam/}_sorted.bam ${file_name};
done;

#for file_name in *sorted.bam.bam;
#do
#mv ${file_name} sorted_bam_files;
#done;

for file_name in *sam.bam;
do
mv ${file_name} ${file_name/sam.bam/}.bam;
done;

for file_name in *_1.fastq;
do
echo mapping ${file_name};
bowtie2 -p 20 -x /dcl01/hansen/data/bjornsson_ht22_chipseq_miseq/ks_chipseq/mm10/mm10 -1 ${file_name} -2 ${file_name/_1.fastq/}_2.fastq -S ${file_name/_1.fastq/}.sam;
done;



for file_name in *sorted.bam;
do
echo removing duplicates from ${file_name};
java -jar /users/lboukas/picard.jar MarkDuplicates I=${file_name} O=${file_name/sorted.bam.bam/}sorted_dedup.bam M=${file_name/_sorted.bam.bam/}.txt REMOVE_DUPLICATES=true;
done;

for file_name in *_sorted.bam;
do
echo removing duplicates from ${file_name};
picard MarkDuplicates I=${file_name} O=${file_name/bam_sorted.bam/}sorted_dedup.bam M=${file_name/bam_sorted.bam/}.txt REMOVE_DUPLICATES=true;
done;

for file_name in *L005.bam_sorted_dedup.bam;
do
echo merging ${file_name/_L005.bam_sorted_dedup.bam};
samtools merge -f -@ 30 ${file_name/_L005.bam_sorted_dedup.bam/}_merged.bam ${file_name/_L005.bam_sorted_dedup.bam/}_L004.bam_sorted_dedup.bam ${file_name} ${file_name/_L005.bam_sorted_dedup.bam/}_L006.bam_sorted_dedup.bam ${file_name/_L005.bam_sorted_dedup.bam/}_L007.bam_sorted_dedup.bam ${file_name/_L005.bam_sorted_dedup.bam/}_L008.bam_sorted_dedup.bam;
done;


for file_name in *L004.sorted_dedup.bam;
do
echo merging ${file_name/_L004.sorted_dedup.bam};
samtools merge -f -@ 30 ${file_name/_L004.sorted_dedup.bam/}_merged.bam ${file_name/_L004.sorted_dedup.bam/}_L001.sorted_dedup.bam ${file_name/_L004.sorted_dedup.bam/}_L002.sorted_dedup.bam ${file_name/_L004.sorted_dedup.bam/}_L003.sorted_dedup.bam;
done;

for file_name in *noMito_merged.bam;
do
echo calling peaks from ${file_name};
macs2 callpeak -t ${file_name} -n ${file_name/_merged.bam/} -g mm --keep-dup all -f BAMPE --bdg;
done;



for file_name in *merged.bam;
do
echo removing mitochondrial reads from ${file_name};
samtools view -@ 30 -b ${file_name} chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY > ${file_name/_merged.bam/}_noMito_merged.bam;
done;

for file_name in *merged.bam;
do
echo indexing ${file_name};
samtools index ${file_name};
done;




###########mapping of rna-seq samples using salmon v10
for file_name in *R1_001.fastq.gz;
do
DYLD_FALLBACK_LIBRARY_PATH=users/lboukas/pub/salmon/salmon-0.10.0_linux_x86_64/lib /users/lboukas/pub/salmon/salmon-0.10.0_linux_x86_64/bin/salmon quant -i /users/lboukas/pub/salmon/salmon-0.10.0_linux_x86_64/bin/mouse_index -l A -1 ${file_name} -2 ${file_name/R1_001.fastq.gz/}R2_001.fastq.gz -p 8 -o quants/${file_name/_R1_001.fastq.gz/};
done;

for file_name in *_1.fastq.gz;
do
DYLD_FALLBACK_LIBRARY_PATH=users/lboukas/pub/salmon/salmon-0.10.0_linux_x86_64/lib /users/lboukas/pub/salmon/salmon-0.10.0_linux_x86_64/bin/salmon quant -i /users/lboukas/human_index -l A -1 ${file_name} -2 ${file_name/_1.fastq.gz/}_2.fastq.gz -p 8 -o quants/${file_name/_1.fastq.gz/};
done;


########fragment length distribution
for file_name in *noMito_merged.bam;
do
echo getting fragment lengths from ${file_name};
samtools view ${file_name} | awk '$9>0' | cut -f 9 | sort | uniq -c | sort -b -k2,2n | sed -e 's/^[ \t]*//' > ${file_name/_merged.bam/}_fragment_length_count.txt;
done;

for file_name in *sizeLessThan150_final.bam;
do
echo indexing ${file_name};
samtools index ${file_name};
done;

########
for file_name in *_noMito_merged.bam;
do
echo getting fragment of defined sizes from ${file_name};
samtools view -@ 30 -h ${file_name} | awk '$9<150 || $1 ~ /^@/' | samtools view -@ 30 -bS - > ${file_name/.bam/}_sizeLessThan150.bam;
done;

for file_name in *sizeLessThan150.bam;
do
echo getting fragment lengths from ${file_name};
samtools view -@ 30 ${file_name} | awk '$9>0' | cut -f 9 | sort | uniq -c | sort -b -k2,2n | sed -e 's/^[ \t]*//' > ${file_name/.bam/}_fragment_length_count.txt;
done;

for file_name in *sizeLessThan150_final.bam;
do
echo calling peaks from ${file_name};
macs2 callpeak -t ${file_name} -n ${file_name/_merged_sizeLessThan150_final.bam/}_onlyLessThan150bpFragments -g mm --keep-dup all -f BAMPE;
done;


for file_name in *all_lessThan150bpFragments_brain_samples_mapq_30.bam;
do
echo calling peaks from ${file_name};
macs2 callpeak --nomodel --nolambda --call-summits -t ${file_name} -n ${file_name/_merged_sizeLessThan150_final.bam/}_onlyLessThan150bpFragments -g mm --keep-dup all -f BAMPE;
done;





#############
for file_name in *noMito_merged.bam;
do
echo downsampling ${file_name};
picard DownsampleSam I=${file_name} O=${file_name/.bam/}_downsampled_25_percent.bam p=0.25 ACCURACY=0.0001 STRATEGY=Chained;
done;

for file_name in *_sizeLessThan150.bam;
do
echo getting fragment of defined sizes from ${file_name};
samtools view -@ 30 -h ${file_name} | awk '$9>-150 || $1 ~ /^@/' | samtools view -@ 30 -bS - > ${file_name/.bam/}_final.bam;
done;



#######
#############
for file_name in *noMito_merged.bam;
do
echo downsampling ${file_name};
samtools view -@ 30 -s 0.25 -b ${file_name} > ${file_name/.bam/}_downsampled_25_percent_with_samtools.bam;
done;





for file_name in *1.fastq;
do
echo mapping ${file_name};
bowtie2 -p 15 -x /dcl01/hansen/data/bjornsson_ht22_chipseq_miseq/ks_chipseq/mm10/mm10 -1 ${file_name} -2 ${file_name/1.fastq/}2.fastq -S ${file_name/1.fastq/}.sam;
done; for file_name in *.sam;
do
samtools view -@ 30 -S -b ${file_name} > ${file_name/.sam/}.bam;
done; for file_name in *.bam;
do
echo sorting ${file_name};
samtools sort -@ 30 -m 2000M ${file_name} ${file_name/.bam/}sorted.bam;
done; for file_name in *sorted.bam.bam;
do
echo removing duplicates from ${file_name};
picard MarkDuplicates I=${file_name} O=${file_name/sorted.bam.bam/}sorted_dedup.bam M=${file_name/sorted.bam.bam/}.txt REMOVE_DUPLICATES=true;
done; for file_name in *dedup.bam;
do
echo getting mitochondrial reads from ${file_name};
samtools view -@ 30 -b ${file_name} chrM > ${file_name/dedup.bam/}_Mito.bam;
done;








for d in */;
do
cd $d; mv *.gz /dcl01/hansen/data/miseq_run_neun_kabuki_03_28_19/all_fastq; cd ..;
done;


for file_name in *sizeLessThan150_final.bam;
do
echo ${file_name}
echo "$(samtools view -c ${file_name})" | bc -l;
done;


for file_name in *_dedup.bam;
do
echo "100 * $(samtools view -c ${file_name} chrM) / $(samtools view -c ${file_name})" | bc -l;
done;

for file_name in *_dedup.bam;
do
echo "100 * $(samtools view -c -f 0x40 -F 0x4 -F 0x8 ${file_name} chrM) / $(samtools view -c -f 0x40 -F 0x4 -F 0x8 ${file_name})" | bc -l;
done;

for file_name in *_dedup.bam;
do
echo "100 * $(samtools view -c -f 0x40 -F 0x4 -F 0x8 ${file_name} chrM) / $(samtools view -c -f 0x40 -F 0x4 -F 0x8 ${file_name})" | bc -l;
done;

for file_name in *_dedup.bam;
do
echo "100 * $(samtools view -c ${file_name} chrM) / $(samtools view -c ${file_name})" | bc -l;
done;

for file_name in *sizeLessThan150_final.bam;
do
echo ${file_name}
echo "$(samtools view -c -f 0x40 -F 0x4 -F 0x8 ${file_name})" | bc -l;
done;


for file_name in *sizeLessThan150_final.bam;
do
echo ${file_name}
echo "$(samtools view -F 0x4 ${file_name})" | cut -f 1 | sort | uniq | wc -l;
done;

for file_name in *_1.fastq;
do
echo mapping ${file_name};
bowtie2 -p 15 -x /dcl01/hansen/data/bjornsson_ht22_chipseq_miseq/ks_chipseq/mm10/mm10 -1 ${file_name} -2 ${file_name/_1.fastq/}_2.fastq -S ${file_name/_1.fastq/}.sam;
done;

for file_name in *dedup.bam;
do
echo removing mitochondrial reads from ${file_name};
samtools view -@ 30 -b ${file_name} chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY > ${file_name/_merged.bam/}_noMito.bam;
done;

for file_name in *dedup.bam;
do
echo indexing ${file_name};
samtools index ${file_name};
done;


#####commands used for merging into a single bam file per condition
#B cells WT
samtools merge -f -@ 30 all_wt_samples_B.bam *B*WT*noMito_merged.bam 21_WT_CD19_ATAC_Cohort_3_S5_noMito_merged.bam 27_CD19_ATAC_KMT2D_Cohort_1_wt_S13_noMito_merged.bam 37_WT_CD19_ATAC_Cohort_3_S4_noMito_merged.bam 26_WT_CD19_ATAC_Cohort_3_S2_noMito_merged.bam 65_CD19_ATAC_KMT2D_Cohort_2_wt_S17_noMito_merged.bam 
#B cells KS1
samtools merge -f -@ 30 all_KS1_samples_B.bam 24_Mut_CD19_ATAC_Cohort_3_S1_noMito_merged.bam 29_Mut_CD19_ATAC_Cohort_3_S3_noMito_merged.bam 31_Mut_CD19_ATAC_Cohort_3_S7_noMito_merged.bam *B*KMT2D*noMito_merged.bam 
#B cells KS2
samtools merge -f -@ 30 all_KS2_samples_B.bam 49_CD19_ATAC_KDM6a_Cohort_1_mut_S6_noMito_merged.bam *B*KDM6a*noMito_merged.bam
#T cells WT
samtools merge -f -@ 30 all_wt_samples_T.bam *T*WT*noMito_merged.bam 27_CD3_ATAC_KMT2D_Cohort_1_wt_S16_noMito_merged.bam 65_CD3_ATAC_KMT2D_Cohort_2_wt_S18_noMito_merged.bam 21_WT_CD90-2_ATAC_Cohort_3_S12_noMito_merged.bam 26_WT_CD90-2_ATAC_Cohort_3_S9_noMito_merged.bam 37_WT_CD90-2_ATAC_Cohort_3_S11_noMito_merged.bam
#T cells KS1
samtools merge -f -@ 30 all_KS1_samples_T.bam 24_Mut_CD90-2_ATAC_Cohort_3_S8_noMito_merged.bam 29_Mut_CD90-2_ATAC_Cohort_3_S10_noMito_merged.bam 31_Mut_CD90-2_ATAC_Cohort_3_S15_noMito_merged.bam *T-KMT2D*noMito_merged.bam
#T cells KS2
samtools merge -f -@ 30 all_KS2_samples_T.bam 49_CD3_ATAC_KDM6a_Cohort_1_mut_S14_noMito_merged.bam *T-KDM6a*noMito_merged.bam

#########RT cohort
#B cells WT 
samtools merge -f -@ 30 all_wt_samples_B_rt_cohort.bam 5-2-6B-ATAC_S2_noMito_merged.bam 5-4--B-ATAC_S4_noMito_merged.bam 5-5-15B-ATAC_S5_noMito_merged.bam 5-6-5B-ATAC_S6_noMito_merged.bam 5-7-44B-ATAC_S7_noMito_merged.bam 6-1-15-2B-ATAC_S9_noMito_merged.bam 6-4-16B-ATAC_S12_noMito_merged.bam
#B cells RT
samtools merge -f -@ 30 all_RT_samples_B_rt_cohort.bam 5-1-42B-ATAC_S1_noMito_merged.bam 5-3-66B-ATAC_S3_noMito_merged.bam 5-8-21B-ATAC_S8_noMito_merged.bam 6-2-8B-ATAC_S10_noMito_merged.bam 6-3-10B-ATAC_S11_noMito_merged.bam
#T cells WT 
samtools merge -f -@ 30 all_wt_samples_T_rt_cohort.bam 6-8--T-ATAC_S16_noMito_merged.bam 7-1-15T-ATAC_S17_noMito_merged.bam 7-2-5T-ATAC_S18_noMito_merged.bam 7-3-44T-ATAC_S19_noMito_merged.bam 7-5-15-2T-ATAC_S21_noMito_merged.bam 7-8-16T-ATAC_S24_noMito_merged.bam 6-6-6T-ATAC_S14_noMito_merged.bam
#T cells RT
samtools merge -f -@ 30 all_RT_samples_T_rt_cohort.bam 6-5-42T-ATAC_S13_noMito_merged.bam 6-7-66T-ATAC_S15_noMito_merged.bam 7-4-21T-ATAC_S20_noMito_merged.bam 7-6-8T-ATAC_S22_noMito_merged.bam 7-7-10T-ATAC_S23_noMito_merged.bam

####
for file_name in *rt_cohort.bam;
do
echo calling peaks from ${file_name};
macs2 callpeak -t ${file_name} -n ${file_name/.bam/} -g mm --keep-dup all -f BAMPE --bdg;
done;

####
for file_name in all*KS2*samples*.bam;
do
echo calling peaks from ${file_name};
macs2 callpeak -t ${file_name} -n ${file_name/.bam/} -g mm --keep-dup all -f BAMPE --bdg;
done;




for file_name in all*KS2*samples*.bam;
do
echo ${file_name};
done;








#####
for file_name in *.fastq;
do
echo mapping ${file_name};
bowtie2 -p 15 -x /dcl01/hansen/data/kmt2d_keratinocytes/chip_seq/hg19/hg19 ${file_name} -S ${file_name/.fastq/}.sam;
done;
#
for file_name in *.fastq;
do
echo mapping ${file_name};
bowtie2 -p 15 -x /dcl01/hansen/data/bjornsson_ht22_chipseq_miseq/ks_chipseq/mm10/mm10 ${file_name} -S ${file_name/.fastq/}.sam;
done;

for file_name in *_1.fastq;
do
echo mapping ${file_name};
bowtie2 -p 15 -x /dcl01/hansen/data/smarca4_study_GSE132293/hg38_index -1 ${file_name} -2 ${file_name/_1.fastq/}_2.fastq -S ${file_name/_1.fastq/}.sam;
done;


samtools merge -f -@ 30 all_wt.bam wt_1.sam.bam_sorted.bamsorted_dedup.bam_noMito.bam wt_2.sam.bam_sorted.bamsorted_dedup.bam_noMito.bam wt_3.sam.bam_sorted.bamsorted_dedup.bam_noMito.bam
samtools merge -f -@ 30 all_mut.bam mut_1.sam.bam_sorted.bamsorted_dedup.bam_noMito.bam mut_2.sam.bam_sorted.bamsorted_dedup.bam_noMito.bam mut_3.sam.bam_sorted.bamsorted_dedup.bam_noMito.bam



###public data atac neurons
samtools merge -f -@ 30 all_samples.bam *noMito.bam

###in house data atac neurons
samtools merge -f -@ 30 all_KS1_samples.bam Cohort8-22*noMito_merged.bam Cohort8-11*noMito_merged.bam Cohort8-9*noMito_merged.bam Cohort8-12*noMito_merged.bam Cohort9-8*noMito_merged.bam Cohort9-14*noMito_merged.bam Cohort9-3*noMito_merged.bam Cohort9-1*noMito_merged.bam Cohort9-2-25*noMito_merged.bam Cohort9-2-20*noMito_merged.bam Cohort9-2-4*noMito_merged.bam Cohort9-2-17*noMito_merged.bam
samtools merge -f -@ 30 all_KS2_samples.bam 10-5-11*noMito_merged.bam 10-5-3*noMito_merged.bam 10-5-9*noMito_merged.bam 10-5-1_*noMito_merged.bam
samtools merge -f -@ 30 all_wt_samples.bam Cohort8-16*noMito_merged.bam Cohort8-19*noMito_merged.bam Cohort8-21*noMito_merged.bam 10-5-2*noMito_merged.bam 10-5-10*noMito_merged.bam 10-5-4*noMito_merged.bam Cohort9-10*noMito_merged.bam Cohort9-6*noMito_merged.bam Cohort9-2_*noMito_merged.bam Cohort9-2-22*noMito_merged.bam Cohort9-2-11*noMito_merged.bam Cohort9-2-12*noMito_merged.bam

macs2 callpeak -t all_KS2_samples.bam -n all_KS2_samples -g mm --keep-dup all -f BAMPE --bdg
macs2 callpeak -t all_KS1_samples.bam -n all_KS1_samples -g mm --keep-dup all -f BAMPE --bdg
macs2 callpeak -t all_wt_samples.bam -n all_wt_samples -g mm --keep-dup all -f BAMPE --bdg

###Jill h3k27me3 chip-seq
samtools merge -f -@ 30 all_ws_samples_chip.bam 705_IP*dedup.bam 719_IP*dedup.bam 720_IP*dedup.bam 729_IP*dedup.bam 739_IP*dedup.bam 748_IP*dedup.bam  
samtools merge -f -@ 30 all_ws_samples_input.bam 705_Input*dedup.bam 719_Input*dedup.bam 720_Input*dedup.bam 729_Input*dedup.bam 739_Input*dedup.bam 748_Input*dedup.bam  

samtools merge -f -@ 30 all_wt_samples_chip.bam 709_IP*dedup.bam 721_IP*dedup.bam 722_IP*dedup.bam 728_IP*dedup.bam 744_IP*dedup.bam 
samtools merge -f -@ 30 all_wt_samples_input.bam 709_Input*dedup.bam 721_Input*dedup.bam 722_Input*dedup.bam 728_Input*dedup.bam 744_Input*dedup.bam 

macs2 callpeak -t 705_IP*dedup.bam -c 705_Input*dedup.bam -n 705 -g mm --keep-dup all -f BAMPE --bdg; 
macs2 callpeak -t 719_IP*dedup.bam -c 719_Input*dedup.bam -n 719 -g mm --keep-dup all -f BAMPE --bdg; 
macs2 callpeak -t 720_IP*dedup.bam -c 720_Input*dedup.bam -n 720 -g mm --keep-dup all -f BAMPE --bdg; 
macs2 callpeak -t 729_IP*dedup.bam -c 729_Input*dedup.bam -n 729 -g mm --keep-dup all -f BAMPE --bdg; 
macs2 callpeak -t 739_IP*dedup.bam -c 739_Input*dedup.bam -n 739 -g mm --keep-dup all -f BAMPE --bdg; 
macs2 callpeak -t 748_IP*dedup.bam -c 748_Input*dedup.bam -n 748 -g mm --keep-dup all -f BAMPE --bdg; 
macs2 callpeak -t 709_IP*dedup.bam -c 709_Input*dedup.bam -n 709 -g mm --keep-dup all -f BAMPE --bdg; 
macs2 callpeak -t 721_IP*dedup.bam -c 721_Input*dedup.bam -n 721 -g mm --keep-dup all -f BAMPE --bdg; 
macs2 callpeak -t 722_IP*dedup.bam -c 722_Input*dedup.bam -n 722 -g mm --keep-dup all -f BAMPE --bdg; 
macs2 callpeak -t 728_IP*dedup.bam -c 728_Input*dedup.bam -n 728 -g mm --keep-dup all -f BAMPE --bdg; 
macs2 callpeak -t 744_IP*dedup.bam -c 744_Input*dedup.bam -n 744 -g mm --keep-dup all -f BAMPE --bdg 

#for drosophila spike in
for file_name in *R1_001.fastq;
do
echo mapping ${file_name};
bowtie2 -p 20 -x /dcl01/hansen/data/fahrner_weaver/weaver_data/chip_seq/spike_in/BDGP6/BDGP6 -1 ${file_name} -2 ${file_name/R1_001.fastq/}R2_001.fastq -S ${file_name/R1_001.fastq/}.sam;
done;


#
MarkDuplicates -I 721_IP_S38_L002.bam_sorted.bam -O 721_IP_S38_L002.bam_sorted.bamsorted_dedup.bam -M 721_IP_S38_L002.bam_sorted.bam.txt -REMOVE_DUPLICATES true




###count number of uniquely mapped reads
for file_name in *dedup.bam;
do
samtools view -F 0x4 ${file_name} | cut -f 1 | sort | uniq | wc -l;
done;

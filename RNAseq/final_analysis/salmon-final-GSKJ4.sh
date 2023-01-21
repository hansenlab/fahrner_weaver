#$ -l mem_free=30G,h_vmem=30G
#$ -l h_fsize=50G
#$ -cwd
#$ -V
#$ -m e
#$ -M cgao13@jhmi.edu

# mapping of rna-seq samples using salmon v1.9.0
{
salmon -v

cd /users/cgao/fahrner_weaver/extdata/RNAseq_GSKJ4

echo "Start time: $(date)"

#echo "Building index"
#salmon index -t gentrome.fa.gz -d decoys.txt -p 8 -i GRCm38_cdna_index_decoy-aware

echo "Mapping reads"
} &>> /users/cgao/fahrner_weaver/RNAseq/final_analysis/log_GSKJ4.txt

for FILE in `ls *.fastq.gz | sed 's/_R[12]_001.fastq.gz//g' | sort -u`
do
	{
	salmon quant -i /users/cgao/fahrner_weaver/extdata/RNAseq/GRCm38_cdna_index_decoy-aware -l A \
		-1 ${FILE}_R1_001.fastq.gz \
		-2 ${FILE}_R2_001.fastq.gz \
		-p 8 --validateMappings --gcBias \
		-o quants_final_GSKJ4/${FILE}_quant
	} 2>> /users/cgao/fahrner_weaver/RNAseq/final_analysis/log_GSKJ4.txt
done

echo "End time: $(date)" &>> /users/cgao/fahrner_weaver/RNAseq/final_analysis/log_GSKJ4.txt

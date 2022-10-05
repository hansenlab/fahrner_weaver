#$ -l mem_free=30G,h_vmem=30G
#$ -l h_fsize=50G
#$ -cwd
#$ -V
#$ -m e
#$ -M cgao13@jhmi.edu

# mapping of rna-seq samples using salmon v1.9.0

salmon -v

cd /users/cgao/fahrner_weaver/extdata/RNAseq

echo "Start time: $(date)"

echo "Building index"
salmon index -t gentrome.fa.gz -d decoys.txt -p 8 -i GRCm38_cdna_index_decoy-aware

echo "Mapping reads"
for FILE in `ls *.fastq.gz | sed 's/_L001_R[12]_001.fastq.gz//g' | sort -u`
do
	salmon quant -i GRCm38_cdna_index_decoy-aware -l A \
		-1 ${FILE}_L001_R1_001.fastq.gz \
		-2 ${FILE}_L001_R2_001.fastq.gz \
		-p 8 --validateMappings --gcBias \
		-o quants_final/${FILE}_quant
done

echo "End time: $(date)"

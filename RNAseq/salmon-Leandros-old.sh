#Christine edited relevant sections from mapping_command_line.R

#$ -l mem_free=5G,h_vmem=5G
#$ -l h_fsize=50G
#$ -cwd
#$ -V
#$ -m e
#$ -M cgao13@jhmi.edu

# mapping of rna-seq samples using salmon v0.10.0 and Leandros' mouse_index

cd /users/cgao/fahrner_weaver/extdata/RNAseq

for file_name in *R1_001.fastq.gz
do
	DYLD_FALLBACK_LIBRARY_PATH=users/lboukas/pub/salmon/salmon-0.10.0_linux_x86_64/lib 
	/users/lboukas/pub/salmon/salmon-0.10.0_linux_x86_64/bin/salmon quant -i /users/lboukas/pub/salmon/salmon-0.10.0_linux_x86_64/bin/mouse_index -l A \
		-1 ${file_name} \
		-2 ${file_name/R1_001.fastq.gz/}R2_001.fastq.gz \
		-p 8 \
		-o quants_Leandros-old/${file_name/_R1_001.fastq.gz/}
done

# Christine edited relevant sections from mapping_command_line.R

# mapping of rna-seq samples using salmon v10
for file_name in *R1_001.fastq.gz
do
#	DYLD_FALLBACK_LIBRARY_PATH=users/lboukas/pub/salmon/salmon-0.10.0_linux_x86_64/lib 
	/users/lboukas/pub/salmon/salmon-0.10.0_linux_x86_64/bin/salmon quant -i /users/lboukas/pub/salmon/salmon-0.10.0_linux_x86_64/bin/mouse_index -l A \
		-1 ${file_name} \
		-2 ${file_name/R1_001.fastq.gz/}R2_001.fastq.gz \
		-p 8 \
		-o quants/${file_name/_R1_001.fastq.gz/}
done

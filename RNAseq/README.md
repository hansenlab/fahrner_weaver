`quants_Christine`: output from running `salmon-Christine.sh`, log files end in `1322722`. `Salmon v1.9.0`, index Ensembl mouse v102.

`quants_Leandros`: Leandros' output copied from `/dcl01/hansen/data/fahrner_weaver/weaver_data/rna_seq/quants`. `Salmon v0.10.0`, index Ensembl mouse v89. Note these are not from running `salmon-Leandros-old.sh` (which resulted in core dump).

Differential expression analysis was performed with permutations of Leandros' and Christine's quant files and R code. Leandros uses `tximport` while Christine uses `tximeta`.

`RNAseq-merge.R`: reduces lists of differentially expressed genes from these various run permutations into a single dataframe containing the fold changes, `RNAseq-merge-log2FC.csv`.

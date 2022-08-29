library(TxDb.Mmusculus.UCSC.mm10.ensGene)
library(tximport)
library(DESeq2)
library(sva)
library(biomaRt)

library(ggplot2)
library(ggrepel)

library(tximeta)
library(SummarizedExperiment)
library(edgeR)
library(limma)
library(statmod)
library(tidyverse)

# ### BEGIN LEANDROS' IMPORT 
# files <- paste0("~/Desktop/temp/Fahrner lab/Data/Weaver/fahrner_weaver/RNAseq/quants_Leandros/",
#                 list.files("quants_Leandros"), "/quant.sf")
# 
# # files <- paste0("~/Desktop/temp/Fahrner lab/Data/Weaver/fahrner_weaver/RNAseq/quants_Christine/", 
# #                 list.files("quants_Christine"), "/quant.sf")
# 
# names(files) <- c("mut_1-1", "mut_1-2", 
#                   "wt_1-1", "wt_1-2", 
#                   "mut_2-1", "mut_2-2",
#                   "mut_3-1", "mut_3-2",
#                   "wt_2-1", "wt_2-2",
#                   "wt_3-1", "wt_3-2",
#                   "wt_4-1", "wt_4-2",
#                   "mut_4-1", "mut_4-2",
#                   "mut_5-1", "mut_5-2",
#                   "wt_5-1", "wt_5-2",
#                   "mut_6-1", "mut_6-2")
# 
# # creates dataframe w/ tx ID and corresponding gene ID
# txdb <- TxDb.Mmusculus.UCSC.mm10.ensGene
# k <- keys(txdb, keytype = "GENEID")
# df <- AnnotationDbi::select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
# tx2gene <- df[, 2:1] 
# 
# # imports tx counts from quant files and match to genes using tx2gene. 
# txi <- tximport(files, type = "salmon", tx2gene = tx2gene, 
#                 countsFromAbundance = "lengthScaledTPM", ignoreTxVersion = TRUE)
# 
# samples <- read.csv(file='RNAseq-samples.csv')
# sampleTable <- data.frame(condition=factor(samples$Genotype))
# dds <- DESeqDataSetFromMatrix(round(txi$counts), sampleTable, design = ~condition)
# ### END LEANDROS' IMPORT 


### BEGIN CHRISTINE'S IMPORT
samples <- read.csv('RNAseq-samples.csv')

# dir <- '~/Desktop/temp/Fahrner lab/Data/Weaver/fahrner_weaver/RNAseq/quants_Christine'
# files <- file.path(dir, paste(samples$Sample.ID, '_quant', sep=''), 'quant.sf')

dir <- '~/Desktop/temp/Fahrner lab/Data/Weaver/fahrner_weaver/RNAseq/quants_Leandros'
files <- file.path(dir, paste(samples$Sample.ID, '_L001', sep=''), 'quant.sf')

file.exists(files)
coldata <- data.frame(files, names=samples$Sample.ID, condition=samples$Genotype, stringsAsFactors = FALSE)

# import transcript abundances
se <- tximeta(coldata)

# summarize abundances to gene level
gse <- summarizeToGene(se)

# load SummarizedExperiment into DESeqDataSet
dds <- DESeqDataSet(gse, design=~condition)
### END CHRISTINE'S IMPORT


### BEGIN CHRISTINE'S ANALYSIS 
dds$condition <- relevel(dds$condition, ref="WT")
plot <- plotPCA(vst(dds))
plot + geom_label_repel(aes(label=gsub("-.*", "", colnames(dds))), show.legend = FALSE)

# collapse technical replicates
dds$sample <- factor(gsub("-.*", "", colnames(dds)))
dds.coll <- collapseReplicates(dds, dds$sample)
stopifnot(all(rowSums(counts(dds[,which(dds$sample=="705")])) == counts(dds.coll[,1])))
dds.coll.counts <- counts(dds.coll)

# filter genes with low median counts
keep <- rowMedians(counts(dds.coll)) > 10
dds.coll.filtered <- dds.coll[keep,]

#dds.coll.filtered$condition <- relevel(dds.coll.filtered$condition, ref="WT")
plot <- plotPCA(vst(dds.coll.filtered))
plot + geom_label_repel(aes(label=gsub("-.*", "", colnames(dds.coll.filtered))), show.legend = FALSE)

# find surrogate variables for batch effects using sva
dds.coll.filtered.counts <- counts(dds.coll.filtered)
mod <- model.matrix(~condition, colData(dds.coll.filtered))
mod0 <- model.matrix(~1, colData(dds.coll.filtered))
svseq <- svaseq(dds.coll.filtered.counts, mod, mod0)

# append surrogate variables
ddssva <- dds.coll.filtered
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
# ddssva$SV3 <- svseq$sv[,3]
# design(ddssva) <- formula(~SV1 + SV2 + SV3 + condition)
design(ddssva) <- formula(~SV1 + SV2 + condition)

# DESeq2 call
ddssva <- DESeq(ddssva)
dds <- ddssva
save(dds, file='dds_CWG-import.rda')
res <- results(dds)

diffexp.subset <- as.data.frame(res[which(res$padj <0.1),])
### END CHRISTINE'S ANALYSIS


# ### BEGIN LEANDROS' ANALYSIS
# # collapses counts from tech replicates
# counts_mat <- txi$counts
# combined_mat_cols <- lapply(seq(1, 22, by = 2), function(xx)
#   rowSums(counts_mat[, c(xx, xx+1), drop = FALSE]))
# 
# combined_mat <- matrix(unlist(combined_mat_cols), ncol = 11)
# rownames(combined_mat) <- names(combined_mat_cols[[1]])
# colnames(combined_mat) <- gsub("_.*", "", colnames(counts_mat)[seq(1, 22, by = 2)])
# 
# ####run differential analysis
# features_by_samples_mat <- combined_mat
# 
# sample_info <- data.frame(genotype = as.factor(colnames(features_by_samples_mat)))
# 
# sampleTable <- data.frame(condition = factor(sample_info$genotype))
# dds <- DESeqDataSetFromMatrix(round(features_by_samples_mat), sampleTable, design = ~condition)
# 
# # relevel - CWG
# dds$condition <- relevel(dds$condition, ref="wt")
# 
# idx <- rowMedians(counts(dds)) > 10
# dat <- counts(dds)[idx,]
# mod <- model.matrix(~ condition, colData(dds))
# mod0 <- model.matrix(~1, colData(dds))
# svseq <- svaseq(dat, mod, mod0)
# 
# dds <- dds[idx, ]
# ddssva <- dds
# ddssva$SV1 <- svseq$sv[,1]
# ddssva$SV2 <- svseq$sv[,2]
# design(ddssva) <- formula(~ SV1 + SV2 + condition)
# ddssva <- DESeq(ddssva)
# 
# dds <- ddssva
# res <- results(dds)
# diffexp.subset <- as.data.frame(res[which(res$padj <0.1),])

# ### END LEANDROS' ANALYSIS


# gene annotation - CWG
ensembl.id <- rownames(diffexp.subset)

ensembl102.mmusculus <- useEnsembl(biomart='genes', dataset='mmusculus_gene_ensembl', version=102)
mgi <- getBM(attributes=c('ensembl_gene_id','mgi_symbol','mgi_description'),
             filters='ensembl_gene_id',
             values=ensembl.id,
             mart=ensembl102.mmusculus)
diffexp.subset$ensembl_gene_id <- ensembl.id
diffexp.subset <- merge(diffexp.subset, mgi, by='ensembl_gene_id')
diffexp.subset <- diffexp.subset[order(diffexp.subset$log2FoldChange, decreasing=TRUE),]
# write.csv(diffexp.subset, file="Leandros-quant-Leandros-code.csv")
write.csv(diffexp.subset, file="CWG-import.csv")

library(TxDb.Mmusculus.UCSC.mm10.ensGene)
library(tximport)
library(DESeq2)
library(sva)
library(biomaRt)

library(ggplot2)
library(ggrepel)

###ezh2 heterozygous missense variant
# edited file path to run locally - CWG
files <- paste0("~/Desktop/temp/Fahrner lab/Data/Weaver/fahrner_weaver/RNAseq/quants_Leandros/",
                list.files("quants_Leandros"), "/quant.sf")

# files <- paste0("~/Desktop/temp/Fahrner lab/Data/Weaver/fahrner_weaver/RNAseq/quants_Christine/", 
#                 list.files("quants_Christine"), "/quant.sf")

names(files) <- c("mut_1-1", "mut_1-2", 
                  "wt_1-1", "wt_1-2", 
                  "mut_2-1", "mut_2-2",
                  "mut_3-1", "mut_3-2",
                  "wt_2-1", "wt_2-2",
                  "wt_3-1", "wt_3-2",
                  "wt_4-1", "wt_4-2",
                  "mut_4-1", "mut_4-2",
                  "mut_5-1", "mut_5-2",
                  "wt_5-1", "wt_5-2",
                  "mut_6-1", "mut_6-2")

######
# creates dataframe w/ tx ID and corresponding gene ID
txdb <- TxDb.Mmusculus.UCSC.mm10.ensGene
k <- keys(txdb, keytype = "GENEID")
df <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
tx2gene <- df[, 2:1] 

# imports tx counts from quant files and match to genes using tx2gene. 
####### transcripts missing from tx2gene: 26405, 38265
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, 
                countsFromAbundance = "lengthScaledTPM", ignoreTxVersion = TRUE)

samples <- read.csv(file='RNAseq-samples.csv')
sampleTable <- data.frame(condition=factor(samples$Genotype))
dds <- DESeqDataSetFromMatrix(round(txi$counts), sampleTable, design = ~condition)

dds$condition <- relevel(dds$condition, ref="WT")
plot <- plotPCA(vst(dds))
plot + geom_label_repel(aes(label=gsub("-.*", "", colnames(dds))), show.legend = FALSE)

# collapse technical replicates
dds$sample <- factor(gsub("-.*", "", colnames(dds)))
dds.coll <- collapseReplicates(dds, dds$sample)
stopifnot(all(rowSums(counts(dds[,which(dds$sample=="mut_1")])) == counts(dds.coll[,1])))
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
save(dds, file='dds_Frankenstein.rda')
res <- results(dds)

diffexp.subset <- as.data.frame(res[which(res$padj <0.1),])

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
write.csv(diffexp.subset, file="Lq-Frankenstein.csv")

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
# # save(dds, file="dds_Leandros-quant-Leandros-code.rda")
# save(dds, file="dds_Christine-quant-Leandros-code.rda")
# 
# diffexp.subset <- as.data.frame(res[which(res$padj <0.1),])
# 
# # gene annotation - CWG
# ensembl.id <- rownames(diffexp.subset)
# 
# ensembl102.mmusculus <- useEnsembl(biomart='genes', dataset='mmusculus_gene_ensembl', version=102)
# mgi <- getBM(attributes=c('ensembl_gene_id','mgi_symbol','mgi_description'),
#              filters='ensembl_gene_id',
#              values=ensembl.id,
#              mart=ensembl102.mmusculus)
# diffexp.subset$ensembl_gene_id <- ensembl.id
# diffexp.subset <- merge(diffexp.subset, mgi, by='ensembl_gene_id')
# diffexp.subset <- diffexp.subset[order(diffexp.subset$log2FoldChange, decreasing=TRUE),]
# # write.csv(diffexp.subset, file="Leandros-quant-Leandros-code.csv")
# write.csv(diffexp.subset, file="Christine-quant-Leandros-code.csv")

###ezh2 inhibitor
files <- paste0("weaver_data/quants_public/", 
                list.files("weaver_data/quants_public"), "/quant.sf")

names(files) <- c("wt_1", "wt_2", "wt_3", "mut_1", "mut_2", "mut_3")

features_by_samples_mat <- counts_mat
colnames(features_by_samples_mat) <- gsub("_.*", "", colnames(features_by_samples_mat))


####
quartz(file = "weaver_vs_public_1.pdf", height = 2.2, width = 4.4, pointsize = 8, type = "pdf")
par(mfrow =  c(1, 2))
hist(res$pvalue, freq = FALSE, breaks = 40, 
     ylim = c(0, 10), col = alpha("red", 0.5), xlab = "p-values", main = "R684C/+ vs WT", xaxt = 'n', yaxt = 'n', font.main = 1)
axis(2, at = c(0, 4, 8))
axis(1, at = c(0, 0.5, 1))
hist(res2$pvalue, freq = FALSE, breaks = 40, ylim = c(0, 10), 
     col = alpha("cornflowerblue", 0.75), main = "GSK126 vs controls", xlab = "p-values", xaxt = 'n', yaxt = 'n', font.main = 1)
axis(2, at = c(0, 4, 8))
axis(1, at = c(0, 0.5, 1))
dev.off()

quartz(file = "weaver_vs_public_2.pdf", height = 2.2, width = 4.4, pointsize = 8, type = "pdf")
par(mfrow =  c(1, 2))
hist(res2$pvalue[which(rownames(res2) %in% rownames(res)[which(res$padj < 0.1)])], freq = FALSE, breaks = 40, 
     ylim = c(0, 18), col = alpha("red", 0.5), xlab = "p-values", main = "R684C/+ vs WT FDR < 0.1", xaxt = 'n', yaxt = 'n', font.main = 1)
axis(2, at = c(0, 9, 18))
axis(1, at = c(0, 0.5, 1))
hist(res2$pvalue[-which(rownames(res2) %in% rownames(res)[which(res$padj < 0.1)])], freq = FALSE, breaks = 40, ylim = c(0, 18), 
     col = alpha("cornflowerblue", 0.75), main = "R684C/+ vs WT FDR >= 0.1", xlab = "p-values", xaxt = 'n', yaxt = 'n', font.main = 1)
axis(2, at = c(0, 9, 18))
axis(1, at = c(0, 0.5, 1))
dev.off()

###
quartz(file = "weaver_ezh2_targets.pdf", height = 2.2, width = 4.4, pointsize = 8, type = "pdf")
par(mfrow =  c(1, 2))
hist(res_RNA$pvalue[which(rownames(res_RNA) %in% mouse_ezh2_targets_strong)], freq = FALSE, breaks = 40, 
     ylim = c(0, 4.5), col = alpha("red", 0.5), xlab = "p-values", main = "EZH2 target genes", xaxt = 'n', yaxt = 'n', font.main = 1)
axis(2, at = c(0, 2, 4))
axis(1, at = c(0, 0.5, 1))
hist(res_RNA$pvalue[-which(rownames(res_RNA) %in% mouse_ezh2_targets_strong)], freq = FALSE, breaks = 40, ylim = c(0, 4.5), 
     col = alpha("cornflowerblue", 0.75), main = "non-EZH2 target genes", xlab = "p-values", xaxt = 'n', yaxt = 'n', font.main = 1)
axis(2, at = c(0, 2, 4))
axis(1, at = c(0, 0.5, 1))
dev.off()



quartz(file = "weaver_pca_plot_1.pdf", height = 3.5, width = 3.5, pointsize = 8, type = "pdf")
par(mfrow = c(1, 2))
vst_mat <- vst(dds)
plotPCA(vst_mat)
dev.off()
###outlier mutants are mice 705, 719


quartz(file = "weaver_pca_plot_2.pdf", height = 3.5, width = 3.5, pointsize = 8, type = "pdf")
par(mfrow = c(1, 2))
vst_mat <- vst(dds)
plotPCA(vst_mat[, -c(1,3)])
dev.off()

quartz(file = "weaver_pca_plot_3_gsk_treatment.pdf", height = 3.5, width = 3.5, pointsize = 8, type = "pdf")
par(mfrow = c(1, 2))
vst_mat <- vst(dds)
plotPCA(vst_mat)
dev.off()
###
proms_mouse_ezh2 <- proms_mouse[which(proms_mouse$gene_id %in% mouse_ezh2_targets_all)]
proms_mouse_ezh2 <- proms_mouse_ezh2[queryHits(findOverlaps(proms_mouse_ezh2, cpg))]

quartz(file = "ezh2_targets_h3k27me3_peaks.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
plot(1, 4656/10676, pch = 19, col = alpha("red", 0.75), bty = 'l', ylab = "% ezh2 targets w/ H3k27me3 peak",
     main = "", xaxt = 'n', yaxt = 'n', xlim = c(0.8, 2.2), ylim = c(0, 1), xlab = "", cex = 1.25)
points(2, 4994/10672, pch = 19, 
       col = alpha("red", 0.75), cex = 1.25)
abline(v = c(1.5, 2.5), lty = "longdash", col = rgb(0,0,0,0.7))
axis(2, at = c(0.1, 0.50, 0.90), labels = c(10, 50, 90))
axis(1, at = c(1, 2), labels = c("WT", "WS"), cex.axis = 0.8, las = 1)
dev.off()



#####
res_weaver_rnaseq$log2FoldChange <- -res_weaver_rnaseq$log2FoldChange
res_weaver_rnaseq_public$log2FoldChange <- -res_weaver_rnaseq_public$log2FoldChange

res_weaver_rnaseq$gene_symbol <- sapply(rownames(res_weaver_rnaseq), 
                               function(xx) unique(proms_mouse$gene_name[which(proms_mouse$gene_id == xx)]))
res_weaver_rnaseq_public$gene_symbol <- sapply(rownames(res_weaver_rnaseq_public), 
                               function(xx) unique(proms_mouse$gene_name[which(proms_mouse$gene_id == xx)]))


res_weaver_significant <- res_weaver_rnaseq[which(res_weaver_rnaseq$padj < 0.1), ]
res_weaver_significant <- res_weaver_significant[order(res_weaver_significant$pvalue), ]

res_weaver_common <- res_weaver_significant[which(rownames(res_weaver_significant) %in% 
                                                     rownames(res_weaver_rnaseq_public)[which(res_weaver_rnaseq_public$padj < 0.1)]), ]
res_weaver_common$logFC_gsk_inh <- sapply(rownames(res_weaver_common), 
         function(xx) res_weaver_rnaseq_public$log2FoldChange[which(rownames(res_weaver_rnaseq_public) == xx)])

write_csv(as.data.frame(res_weaver_significant), "weaver_rnaseq.csv")
write_csv(as.data.frame(res_weaver_common), "weaver_rnaseq_shared_with_gsk_inh.csv")


###
res_weaver_rna_and_chip <- qvalue(res_weaver_rnaseq$pvalue[which(rownames(res_weaver_rnaseq) %in% 
                                                            res_weaver_significant_chipseq_proms_df$gene_id)], pi0.method = "bootstrap", 
                                  fdr.level = 0.1)
res_weaver_rna_and_chip <- res_weaver_rnaseq[which(rownames(res_weaver_rnaseq) %in% 
                                                     res_weaver_significant_chipseq_proms_df$gene_id)[
                                                       which(res_weaver_rna_and_chip$significant == TRUE)], ]

res_weaver_rna_and_chip$logFC_chip <- sapply(rownames(res_weaver_rna_and_chip), function(xx) 
  res_proms_weaver$log2FoldChange[which(res_proms_weaver$gene_id == xx)])

res_weaver_rna_and_chip$mean_logFC_chip <- sapply(1:dim(res_weaver_rna_and_chip)[1], function(xx) 
  mean(unlist(res_weaver_rna_and_chip$logFC_chip[xx])))

res_weaver_rna_and_chip$pval_chip <- sapply(rownames(res_weaver_rna_and_chip), function(xx) 
  res_proms_weaver$pvalue[which(res_proms_weaver$gene_id == xx)])


res_weaver_rna_and_chip$mean_pval_chip <- sapply(1:dim(res_weaver_rna_and_chip)[1], function(xx) 
  mean(unlist(res_weaver_rna_and_chip$pval_chip[xx])))

quartz(file = "weaver_chip_logFC.pdf", height = 2.2, width = 6.5, pointsize = 8, type = "pdf")
par(mfrow = c(1, 3))
hist(res_granges_weaver$log2FoldChange, 
     freq = FALSE, breaks = 40, 
     col = alpha("red", 0.57), xlab = "log2FC", main = "ChIP-seq\n(all peaks)", 
     font.main = 1, lty = 0, xaxt = 'n', yaxt = 'n')
abline(v = 0, lty = "longdash", lwd = 2.5)
axis(2, at = c(0, 1.5))
axis(1, at = c(-1, 0, 1))
hist(res_granges_weaver$log2FoldChange[which(res_granges_weaver$padj < 0.1)], 
     freq = FALSE, breaks = 40, 
     ylim = c(0, 4.5), col = alpha("red", 0.57), xlab = "log2FC", main = "ChIP-seq\n(significant peaks)", 
     xaxt = 'n', yaxt = 'n', font.main = 1, lty = 0, xlim = c(0, 1))
axis(2, at = c(0, 4))
axis(1, at = c(0, 0.5, 1))
hist(res_proms_weaver$log2FoldChange[which(res_proms_weaver$padj < 0.1)], 
     freq = FALSE, breaks = 35, 
     ylim = c(0, 5), col = alpha("red", 0.57), xlab = "log2FC", main = "ChIP-seq\n(significant peaks)", 
     xaxt = 'n', yaxt = 'n', font.main = 1, lty = 0, xlim = c(0, 1))
axis(2, at = c(0, 4))
axis(1, at = c(0, 0.5, 1))
dev.off()

res_weaver_rna_and_chip$gene_id <- rownames(res_weaver_rna_and_chip)
write_csv(as.data.frame(res_weaver_rna_and_chip), "weaver_genes_with_expression_and_h3k27me3_disruption.csv")


#####
quartz(file = "weaver_rnaseq_volcano_plot.pdf", height = 4.4, width = 6.5, pointsize = 8, type = "pdf")
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = c(''),
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = max(res$pvalue[which(res$padj < 0.05)]),
                FCcutoff = 1,
                pointSize = 1.5,
                labSize = 6.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                colAlpha = 4/5,
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0, colConnectors = 'black')

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = c(''),
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = max(res$pvalue[which(res$padj < 0.05)]),
                FCcutoff = 1,
                hlineType = "blank",
                pointSize = 1.5,
                labSize = 6.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                colAlpha = 4/5,
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0, colConnectors = 'black')
dev.off()




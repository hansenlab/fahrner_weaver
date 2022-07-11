load(file = "weaver_data/count_reads_in_combined_peaks_from_condition_specific_bam_files_neurons_weaver.rda")

peaks_by_samples <- count_reads_in_peaks$counts
sample_info <- data.frame(genotype = rep("WT", 11), stringsAsFactors = FALSE)
sample_info$genotype[c(1, 3, 4, 8, 9, 11)] <- "WS"

sampleTable <- data.frame(condition = factor(sample_info$genotype))
dds <- DESeqDataSetFromMatrix(round(peaks_by_samples), sampleTable, design = ~condition)

#compute size factors by calculating the number of reads mapping to drosophila genome after excluding duplicated reads
#I do this only for the IP samples with the following command:
#samtools view -F 0x4 ${file_name} | cut -f 1 | sort | uniq | wc -l;

mapped_reads <- c(7503, 8230, 10835, 10410, 8756, 8029, 7585, 9301, 8584, 7717, 10141) 
norm_factors <- min(mapped_reads)/mapped_reads
sizeFactors(dds) <- norm_factors

#sampleTable <- data.frame(condition = factor(sample_info$genotype[-4])) #because sample 4 is an outlier
#dds <- DESeqDataSetFromMatrix(round(peaks_by_samples[, -4]), sampleTable, design = ~condition)

idx <- rowMedians(counts(dds)) > 10 & rowMedians(counts(dds)) < 5000 
dat <- counts(dds)[idx,]
mod <- model.matrix(~ condition, colData(dds))
mod0 <- model.matrix(~1, colData(dds))
svseq <- svaseq(dat, mod, mod0)

dds <- dds[idx, ]
ddssva <- dds
ddssva$SV1 <- svseq$sv[,1]
design(ddssva) <- formula(~ SV1 + condition)
ddssva <- DESeq(ddssva)


dds <- ddssva
res <- results(dds)
#if (length(which(is.na(res$padj))) > 0){res <- res[-which(is.na(res$padj)), ]} ###too many adjusted p values are NAs in KS1. investigate this!!!
res$P.Value <- res$pvalue



###
peak_df <- count_reads_in_peaks$annotation
rownames(peak_df) <- peak_df$GeneID

res$log2FoldChange <- -res$log2FoldChange
res$chr <- peak_df[rownames(res), "Chr"]
res$start <- peak_df[rownames(res), "Start"]
res$end <- peak_df[rownames(res), "End"]
res_granges_weaver <- makeGRangesFromDataFrame(res, keep.extra.columns = TRUE)
prom_overlaps <- findOverlaps(res_granges_weaver, proms_mouse)
res_proms <- res_granges_weaver[queryHits(prom_overlaps)]
res_proms$gene_id <- proms_mouse$gene_id[subjectHits(prom_overlaps)]
res_proms$gene_name <- proms_mouse$gene_name[subjectHits(prom_overlaps)]
res_proms_weaver <- res_proms[-which(duplicated(res_proms))]

peaks_by_prom <- split(res_proms, res_proms$gene_id)

res_proms_diff <- res_proms[which(res_proms$padj < 0.1)]
diff_peaks_by_prom <- split(res_proms_diff, res_proms_diff$gene_id)

logfc_chip <- sapply(diff_peaks_by_prom, function(xx) mean(xx$log2FoldChange))
res_RNA_df <- as.data.frame(res_RNA)
res_RNA_diff <- res_RNA_df[which(res_RNA_df$padj < 0.1), ]

logfc_chip <- logfc_chip[rownames(res_RNA_diff)]
logfc_rna <- res_RNA_df[names(logfc_chip), "log2FoldChange"]


###

quartz(file = "h3k27me3_vs_expression_logfc.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
plot(-logfc_chip, -logfc_rna, 
     cex = 1.15, col = alpha("red", 0.4), xlab = "promoter H3K27me3 logFC", ylab = "gene expression logFC", 
     main = "", font.main = 1, pch = 19, xlim = c(0, 1), ylim = c(-2.5, 2), bty = 'l', xaxt = 'n', yaxt = 'n', cex.lab = 1.25)
axis(1, at = c(0, 0.5, 1), cex.axis = 1.2)
axis(2, at = c(-2, 0, 2), cex.axis = 1.2)
abline(h=0, lty = "longdash", col = rgb(0,0,0,0.4))
dev.off()


quartz(file = "h3k27me3_vs_expression_pval.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
plot(density(res_RNA$pvalue[which(rownames(res_RNA) %in% 
             res_proms$gene_id[which(res_proms$pvalue <= quantile(res_proms$pvalue, 0.1))])], from = 0, to = 1, bw = 0.035), 
     ylim = c(0, 5), col = alpha("red", 0.57), xlab = "p-values", lwd = 2.5, xaxt = 'n', yaxt = 'n', main = "")
lines(density(res_RNA$pvalue[-which(rownames(res_RNA) %in% 
             res_proms$gene_id[which(res_proms$pvalue <= quantile(res_proms$pvalue, 0.1))])], from = 0, to = 1), 
      col = "cornflowerblue", lwd = 2.5)
axis(2, at = c(0, 2.5, 5))
axis(1, at = c(0, 0.5, 1))
dev.off()


#####
wt_samples <- getGrangesFromNarrowPeaks('weaver_data/all_wt_samples_peaks.narrowPeak')
ws_samples <- getGrangesFromNarrowPeaks('weaver_data/all_ws_samples_peaks.narrowPeak')

overlaps <- findOverlaps(wt_samples, ws_samples)

wt_samples <- wt_samples[queryHits(overlaps)]
ws_samples <- ws_samples[subjectHits(overlaps)]

###
median_p <- sapply(seq(1000, 5000, by = 250), function(xx) 
  median(res_weaver_rna$pvalue[which(rownames(res_weaver_rna) %in% unlist(res_proms$gene_id[which(rank(res_proms$pvalue) < xx)]))]))


###
quartz(file = "weaver_pca_plot_chip_1.pdf", height = 3.5, width = 3.5, pointsize = 8, type = "pdf")
#par(mfrow = c(1, 2))
vst_mat <- vst(dds)
plotPCA(vst_mat)
dev.off()

quartz(file = "weaver_pca_plot_chip_2.pdf", height = 3.5, width = 3.5, pointsize = 8, type = "pdf")
#par(mfrow = c(1, 2))
vst_mat <- vst(dds)
plotPCA(vst_mat)
dev.off()



###
quartz(file = "weaver_rna_vs_chip_pval_histogram.pdf", height = 2.2, width = 4.4, pointsize = 8, type = "pdf")
par(mfrow =  c(1, 2))
hist(res_weaver_rna$pvalue, freq = FALSE, breaks = 40, 
     ylim = c(0, 3), col = alpha("red", 0.5), xlab = "p-values", main = "rna-seq", xaxt = 'n', yaxt = 'n', font.main = 1)
axis(2, at = c(0, 1.5, 3))
axis(1, at = c(0, 0.5, 1))
hist(res$pvalue, freq = FALSE, breaks = 40, ylim = c(0, 1.2), 
     col = alpha("cornflowerblue", 0.75), main = "chip-seq (with outlier)", xlab = "p-values", xaxt = 'n', yaxt = 'n', font.main = 1)
axis(2, at = c(0, 0.5, 1))
axis(1, at = c(0, 0.5, 1))
dev.off()

quartz(file = "weaver_rna_vs_chip_pval_histogram_2.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
hist(res_no_outlier$pvalue, freq = FALSE, breaks = 40, 
     ylim = c(0, 1.2), col = alpha("cornflowerblue", 0.5), xlab = "p-values", main = "chip-seq (no ourlier)", xaxt = 'n', yaxt = 'n', font.main = 1)
axis(2, at = c(0, 0.5, 1))
axis(1, at = c(0, 0.5, 1))

dev.off()


###
quartz(file = "weaver_ezh2_targets_chip.pdf", height = 2.2, width = 4.4, pointsize = 8, type = "pdf")
par(mfrow =  c(1, 2))
hist(res_proms$pvalue[which(res_proms$gene_id %in% mouse_ezh2_targets_strong)], freq = FALSE, breaks = 40, 
     ylim = c(0, 1.5), col = alpha("red", 0.5), xlab = "p-values", main = "EZH2 target genes", xaxt = 'n', yaxt = 'n', font.main = 1)
axis(2, at = c(0, 0.6, 1.2))
axis(1, at = c(0, 0.5, 1))
hist(res_proms$pvalue[-which(res_proms$gene_id %in% mouse_ezh2_targets_strong)], freq = FALSE, breaks = 40, ylim = c(0, 4.5), 
     col = alpha("cornflowerblue", 0.75), main = "non-EZH2 target genes", xlab = "p-values", xaxt = 'n', yaxt = 'n', font.main = 1)
axis(2, at = c(0, 1.5))
axis(1, at = c(0, 0.6, 1.2))
dev.off()


###
quartz(file = "weaver_chip_pval_histogram_spike_in_normalization.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
hist(res$pvalue, freq = FALSE, breaks = 40, 
     ylim = c(0, 3.2), col = alpha("red", 0.5), xlab = "p-values", main = "chip-seq (spike-in normalization)", 
     xaxt = 'n', yaxt = 'n', font.main = 1)
axis(2, at = c(0, 1.5, 3))
axis(1, at = c(0, 0.5, 1))
dev.off()

quartz(file = "h3k27me3_vs_expression_pval.pdf.pdf", height = 2.2, width = 4.4, pointsize = 8, type = "pdf")
par(mfrow = c(1, 2))
hist(res_RNA$pvalue[which(rownames(res_RNA) %in% res_proms$gene_id[which(res_proms$pvalue <= quantile(res_proms$pvalue, 0.1))])], 
     freq = FALSE, breaks = 40, 
     ylim = c(0, 5), col = alpha("red", 0.57), xlab = "p-values", main = "RNA-seq\n(promoters w/ differential\nchip signal)", 
     xaxt = 'n', yaxt = 'n', font.main = 1, lty = 0)
axis(2, at = c(0, 2.5, 5))
axis(1, at = c(0, 0.5, 1))
hist(res_RNA$pvalue[-which(rownames(res_RNA) %in% res_proms$gene_id[which(res_proms$pvalue <= quantile(res_proms$pvalue, 0.1))])], 
     freq = FALSE, breaks = 40, 
     ylim = c(0, 5), col = "cornflowerblue", xlab = "p-values", main = "RNA-seq\n(promoters w/out differential\nchip signal)", 
     xaxt = 'n', yaxt = 'n', font.main = 1, lty = 0)
axis(2, at = c(0, 2.5, 5))
axis(1, at = c(0, 0.5, 1))
dev.off()


###check individual samples to see if there's anything obviously wrong with the outlier
ranges_720 <- getGrangesFromNarrowPeaks('weaver_data/720_peaks.narrowPeak')
test <- findOverlaps(proms_mouse[which(proms_mouse$gene_id %in% mouse_ezh2_targets_strong)], ranges_720)
length(unique(queryHits(test)))

#####
library(Repitools)
res_weaver_significant_chipseq_all <- res_granges_weaver[which(res_granges_weaver$padj < 0.1), ]
res_weaver_significant_chipseq_proms <- res_proms_weaver[which(res_proms_weaver$padj < 0.1), ]

res_weaver_significant_chipseq_all_df <- annoGR2DF(res_weaver_significant_chipseq_all)
res_weaver_significant_chipseq_proms_df <- annoGR2DF(res_weaver_significant_chipseq_proms)
res_weaver_significant_chipseq_non_proms_df <- res_weaver_significant_chipseq_all_df[-which(
  rownames(res_weaver_significant_chipseq_all_df) %in% rownames(res_weaver_significant_chipseq_proms_df)), ]


write_csv(res_weaver_significant_chipseq_non_proms_df, "weaver_chipseq_outside_promoters.csv")
write_csv(res_weaver_significant_chipseq_proms_df, "weaver_chipseq_promoters.csv")



###
res_weaver_significant_chipseq_non_proms_grange <- makeGRangesFromDataFrame(res_weaver_significant_chipseq_non_proms_df, 
                                                                            keep.extra.columns = TRUE)

res_weaver_significant_chipseq_non_proms_grange <- res_weaver_significant_chipseq_non_proms_grange + 50000
prom_overlaps2 <- findOverlaps(res_weaver_significant_chipseq_non_proms_grange, proms_mouse)
res_proms2 <- res_weaver_significant_chipseq_non_proms_grange[queryHits(prom_overlaps2)]
res_proms2$gene_id <- proms_mouse$gene_id[subjectHits(prom_overlaps2)]
res_proms2$gene_name <- proms_mouse$gene_name[subjectHits(prom_overlaps2)]
res_proms_weaver_2 <- res_proms2[-which(duplicated(res_proms2))]


#
quartz(file = "rna_vs_distant_chip.pdf", height = 2.2, width = 4.4, pointsize = 8, type = "pdf")
par(mfrow =  c(1, 2))
hist(res_weaver_rnaseq$pvalue[which(rownames(res_weaver_rnaseq) %in% unlist(res_proms_weaver_2$gene_id))], freq = FALSE, breaks = 40, 
     ylim = c(0, 4), col = alpha("red", 0.5), xlab = "p-values", main = "overlap with ChIP-seq\nsignificant peaks", xaxt = 'n', yaxt = 'n', font.main = 1)
axis(2, at = c(0, 2, 4))
axis(1, at = c(0, 0.5, 1))
hist(res_weaver_rnaseq$pvalue[-which(rownames(res_weaver_rnaseq) %in% unlist(res_proms_weaver_2$gene_id))], freq = FALSE, breaks = 40, ylim = c(0, 4), 
     col = alpha("cornflowerblue", 0.75), main = "no overlap with ChIP-seq\nsignificant peaks", xlab = "p-values", xaxt = 'n', yaxt = 'n', font.main = 1)
axis(2, at = c(0, 2, 4))
axis(1, at = c(0, 0.5, 1))
dev.off()







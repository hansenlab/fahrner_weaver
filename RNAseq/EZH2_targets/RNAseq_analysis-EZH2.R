library(rtracklayer)
library(EnsDb.Hsapiens.v75)
library(BSgenome.Hsapiens.UCSC.hg19)
library(biomaRt)
library(DESeq2)

# session <- browserSession('UCSC')
# genome(session) <- 'hg19'
# # query <- ucscTableQuery(session, 'Txn Factor ChIP', )
# # tableName(query) <- 'wgEncodeRegTfbsClusteredV3'
# query <- ucscTableQuery(session, table='wgEncodeRegTfbsClusteredV3')
# tf_table_query <- getTable(query)
# tf_table_query_EZH2 <- tf_table_query[which(tf_table_query$name=='EZH2'),]

bed <- as.data.frame(read.table("wgEncodeRegTfbsClusteredWithCellsV3.bed",
                                header = FALSE, sep="\t", stringsAsFactors=FALSE, quote=""))
colnames(bed) <- c('chrom','chromStart','chromEnd','factor','score','cell')

bed$sourceCount <- lengths(str_split(bed$cell, ','))
bed_EZH2 <- bed[which(bed$factor=='EZH2'),]
table(unlist(str_split(bed_EZH2$cell, ',')))

# get GRanges of human EZH2 binding peaks
tf_table <- read.csv(file='tf_table.csv')
tf_table_EZH2 <- tf_table[which(tf_table$name=='EZH2'),]
EZH2_h19 <- tf_table_EZH2[,c('chrom','chromStart','chromEnd','name','sourceCount')]
EZH2_h19_granges <- makeGRangesFromDataFrame(EZH2_h19, keep.extra.columns=TRUE)
genome(EZH2_h19_granges) <- 'hg19'
save(EZH2_h19_granges, file='EZH2_h19_granges.rda')

# separate GRanges by level of evidence
load(file='EZH2_h19_granges.rda')
EZH2_all <- EZH2_h19_granges
EZH2_mod <- EZH2_h19_granges[which(EZH2_h19_granges$sourceCount >= 2)]
EZH2_strong <- EZH2_h19_granges[which(EZH2_h19_granges$sourceCount >= 5)]

# get GRanges of all human promoters
edb <- EnsDb.Hsapiens.v75
proms_all <- promoters(edb, upstream=2000, downstream=2000,
                       columns = c("gene_name", "tx_id", "tx_cds_seq_start", "tx_cds_seq_end",
                                   "tx_biotype", "gene_id"))
genome(seqinfo(proms_all)) <- 'hg19'
seqlevelsStyle(proms_all) <- 'ucsc'

# overlap human EZH2 binding peaks with human promoters (human EZH2 targets)
proms_all <- proms_all[which(seqnames(proms_all) %in% names(Hsapiens)[1:24])]
EZH2_target_id_all <- unique(proms_all$gene_id[queryHits(findOverlaps(proms_all, EZH2_all))])
EZH2_target_id_mod <- unique(proms_all$gene_id[queryHits(findOverlaps(proms_all, EZH2_mod))])
EZH2_target_id_strong <- unique(proms_all$gene_id[queryHits(findOverlaps(proms_all, EZH2_strong))])

# find mouse homologs of human EZH2 targets
human <- useMart('ensembl', dataset='hsapiens_gene_ensembl')
searchAttributes(human, pattern='mmusculus')
human_mouse_homologs <- getBM(attributes = c("ensembl_gene_id", "external_gene_name",
                                             "mmusculus_homolog_ensembl_gene", 
                                             "mmusculus_homolog_associated_gene_name", 
                                             "mmusculus_homolog_orthology_confidence", 
                                             "mmusculus_homolog_perc_id_r1"),
                              filters = "ensembl_gene_id",
                              values = EZH2_target_id_all,
                              mart = human)
human_mouse_homologs <- human_mouse_homologs[
  which(human_mouse_homologs$mmusculus_homolog_orthology_confidence==1),]
mouse_EZH2_target_all <- unique(human_mouse_homologs$mmusculus_homolog_ensembl_gene)
save(mouse_EZH2_target_all, file='mouse_EZH2_target_all.rda')

human_mouse_homologs_strong <- getBM(attributes = c("ensembl_gene_id", "external_gene_name",
                                             "mmusculus_homolog_ensembl_gene", 
                                             "mmusculus_homolog_associated_gene_name", 
                                             "mmusculus_homolog_orthology_confidence", 
                                             "mmusculus_homolog_perc_id_r1"),
                              filters = "ensembl_gene_id",
                              values = EZH2_target_id_strong,
                              mart = human)
human_mouse_homologs_strong <- human_mouse_homologs_strong[
  which(human_mouse_homologs_strong$mmusculus_homolog_orthology_confidence==1),]
mouse_EZH2_target_strong <- unique(human_mouse_homologs_strong$mmusculus_homolog_ensembl_gene)
save(mouse_EZH2_target_strong, file='mouse_EZH2_target_strong.rda')

# p-value histograms of EZH2 targets in final_analysis data
setwd("~/Desktop/temp/Fahrner lab/Data/Weaver/fahrner_weaver/RNAseq/final_analysis/BigStuff")
load(file='dds_final.rda')
res <- results(dds)
setwd("~/Desktop/temp/Fahrner lab/Data/Weaver/fahrner_weaver/RNAseq/EZH2_targets")
load(file='mouse_EZH2_target_all.rda')
load(file='mouse_EZH2_target_strong.rda')

setwd("~/Desktop/temp/Fahrner lab/Data/Weaver/fahrner_weaver/RNAseq/EZH2_targets/plots")
setEPS()
postscript(file='221120_p-val-EZH2.eps')
hist(res$pvalue[which(rownames(res) %in% mouse_EZH2_target_all)], 
     freq=FALSE, breaks=50,
     xlab='p-values (R684C/+ vs WT)',
     main='EZH2 target genes',
     ylim=c(0,4.5),
     col='cornflowerblue',
     cex.lab=1.5,
     cex.axis=1.5,
     cex.main=1.75)
dev.off()

setEPS()
postscript(file='221120_p-val-non-EZH2.eps')
hist(res$pvalue[!(rownames(res) %in% mouse_EZH2_target_all)], 
     freq=FALSE, breaks=50,
     xlab='p-values (R684C/+ vs WT)',
     main='non-EZH2 target genes',
     ylim=c(0,4.5),
     col='gray60',
     cex.lab=1.5,
     cex.axis=1.5,
     cex.main=1.75)
dev.off()

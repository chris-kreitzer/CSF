##----------------+
## Cnlr against baseline
##----------------+
## start: 08/31/2022
## revision: 10/20/2022
## chris-kreitzer


clean()
gc()
.rs.restartR()
source('~/Documents/GitHub/CSF/Scripts/UtilityFunctions.R')
setup(working.path = '~/Documents/MSKCC/Subhi/CSF/C_006884/')
library(patchwork)
library(data.table)
library(dplyr)

# Test on cnlr against baseline
segs = facets_out$segs
dipLogR = facets_out$dipLogR
snps = facets_out$snps

cn0_dipLogR = unique(segs$cnlr.median.clust)[which.min(abs(unique(segs$cnlr.median.clust)-dipLogR))]
cn0_segs = segs[segs$cnlr.median.clust == cn0_dipLogR, ]
cn0_snps = snps[snps$segclust %in% cn0_segs$segclust, ]
cn0_snps = cn0_snps[between(cn0_snps$cnlr, quantile(cn0_snps$cnlr, .25), quantile(cn0_snps$cnlr, .75)), ] # remove noise

cnlr = snps[which(snps$chrom == 7 & snps$maploc >= 55086714 & snps$maploc <= 55324313), ]
ifelse(mean(cnlr$cnlr) > mean(cn0_snps$cnlr),
       two_sample_z(cn0_snps$cnlr, cnlr$cnlr),
       two_sample_z(cnlr$cnlr, cn0_snps$cnlr))
fold_change = ifelse(mean(cnlr$cnlr) - mean(cn0_snps$cnlr) < 0,
                     -2^(-(mean(cnlr$cnlr) - mean(cn0_snps$cnlr))),
                     2^(mean(cnlr$cnlr) - mean(cn0_snps$cnlr)))


# Perform Z test, calculate fold change
genes_all[, mean_cnlr := mean(unlist(cnlr)), by = seq_len(nrow(genes_all))][, `:=` (
  pval = ifelse(mean_cnlr > mean(cn0_snps$cnlr),
                two_sample_z(cn0_snps$cnlr, unlist(cnlr)),
                two_sample_z(unlist(cnlr), cn0_snps$cnlr)),
  fold_change = ifelse(mean_cnlr - mean(cn0_snps$cnlr) < 0,
                       -2^(-(mean_cnlr - mean(cn0_snps$cnlr))),
                       2^(mean_cnlr - mean(cn0_snps$cnlr)))
  ), by = seq_len(nrow(genes_all))]


two_sample_z = function(a, b) {
    if (length(a) < 5 | length(b) < 5) {
        NA_real_
    } else {
        se_a = sd(a)/sqrt(length(a))
        se_b = sd(b)/sqrt(length(b))
        se = sqrt(se_a^2 + se_b^2)
        z = (mean(a)-mean(b))/se
        pnorm(z, lower.tail = TRUE)
    }
}






##----------------+
## EGFR amplification 
## based on FACETS:
##----------------+

setwd('~/Documents/MSKCC/Subhi/CSF/')

binary = read.csv('Data/FINAL_samples/CSF_binary_QC_true.txt', sep = '\t')
View(binary[row.names(binary), 'P.0033342.T01.IM6'])
cbind(row.names(binary), binary$P.0033342.T01.IM6)

EGFR = binary[which(row.names(binary) == 'EGFR'), ]
View(EGFR)









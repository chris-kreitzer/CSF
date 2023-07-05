## Cell Lines Subhi
## 07/05/2023

setwd('~/Documents/MSKCC/11_CSF/')


library(patchwork)
library(dplyr)

source('~/Documents/GitHub/CSF/Scripts/UtilityFunctions.R')
source('~/Documents/GitHub/CSF/Scripts/FacetsPlot.R')
source('~/Documents/GitHub/CSF/Scripts/cf_Plot.R')
source('~/Documents/GitHub/CSF/Scripts/ReadSnpMatrix.R')
source('~/Documents/GitHub/CSF/Scripts/gene_closeup.R')
source('~/Documents/GitHub/CSF/Scripts/GMM.R')
source('~/Documents/GitHub/CSF/Scripts/ExonicMapping_CDKN2A.R')


GOIs = c('CDKN2A','CDK4','CDK6','PTEN','EGFR','PDGFRA',
         'KIT','KDR','MET','MDM2','MDM4','RB1','NF1',
         'TP53','FGF4', 'FGF19', 'PIK3CA', 'BRAF')


path = '09_Cell_Line/FROZENPOOLEDNORMAL_IMPACT505_V2.rg.md.abra.printreads__s_TS543_14452/FROZENPOOLEDNORMAL_IMPACT505_V2.rg.md.abra.printreads__s_TS543_14452.rg.md.abra.printreads.dat.gz'
sample = 'FROZENPOOLEDNORMAL_IMPACT505_V2.rg.md.abra.printreads__s_TS543_14452'

countmatrix = readsnpmatrix(path = path)
countmatrix = as.data.frame(countmatrix)

snps = facetsSuite::run_facets(read_counts = countmatrix,
                               cval = 150,
                               dipLogR = NULL,
                               snp_nbhd = 250,
                               seed = 100, 
                               genome = 'hg19', 
                               ndepth = 25)
het_snps = snps$snps[which(snps$snps$het == 1), ]

dir.create(path = paste0('09_Cell_Line/', sample))
pdf(file = paste0('09_Cell_Line/', sample, '/', 'Normal_Het_Distribution.pdf'), width = 4.5, height = 4.5)
plot(density(het_snps$vafN[which(het_snps$het == 1)]), main = 'allele freq. het. SNPs\nmatched NORMAL')
dev.off()
norm_density = density(het_snps$vafN)
norm_density_max = norm_density$x[which.max(norm_density$y)]
ifelse(norm_density_max <= 0.4 | norm_density_max >= 0.6, print('STOP'), print('OK'))
pdf(file = paste0('09_Cell_Line//', sample, '/', 'TUMOR_Het_Distribution.pdf'), width = 4.5, height = 4.5)
plot(density(het_snps$vafT[which(het_snps$het == 1)]), main = 'allele freq. het. SNPs\nTUMOR')
dev.off()
rm(snps, het_snps)
sample

##--
fourth = facetsSuite::run_facets(read_counts = countmatrix,
                                 cval = 100,
                                 dipLogR = NULL,
                                 snp_nbhd = 50,
                                 seed = 100, 
                                 genome = 'hg19', 
                                 ndepth = 25)

snps = fourth$snps

##-- run reduced fits (alias 'imitated' MSK-pipeline (me))
##-- Assume that all the CnLR are normally distributed
Mean = mean(snps$cnlr, na.rm = T)
Sd = sd(snps$cnlr, na.rm = T)

x = seq(floor(min(snps$cnlr)), ceiling(max(snps$cnlr)), 0.01)
xp = pnorm(x, mean = Mean, sd = Sd)

##-- lower (10%) and upper 90% threshold
lower = x[max(which(xp <= 0.1, arr.ind = T))]
upper = x[min(which(xp >= 0.9, arr.ind = T))]

##-- run per gene of interest (GOI)
chris = data.frame()
for(gene in unique(GOIs)){
  position = snps[which(snps$chrom == genes_hg19$chrom[which(genes_hg19$gene == gene)] &
                          snps$maploc >= genes_hg19$start[which(genes_hg19$gene == gene)] &
                          snps$maploc <= genes_hg19$end[which(genes_hg19$gene == gene)]), ]
  gene_cnlr_mean = mean(position$cnlr)
  gene_call = ifelse(gene_cnlr_mean >= upper, 2, 
                     ifelse(gene_cnlr_mean <= lower, -2, 0))
  hugo = gene
  
  gene_out = data.frame(id = sample,
                        gene = hugo,
                        gene_cnlr_mean = gene_cnlr_mean,
                        call = gene_call,
                        upper = upper,
                        lower = lower,
                        n_hets = length(snps$maploc[which(snps$het == 1)]))
  
  chris = rbind(chris, gene_out)
}

chris

pdf(file = paste0('09_Cell_Line//', sample, '/', 'IMPACT-like.CnLR_distribution.pdf'), 
    width = 6, height = 6)
plot(x, 
     pnorm(x, mean = Mean, sd = Sd), type = "l",
     ylim = c(0, 1), ylab = "P(X < x)", lwd = 2, col = "red")
##-- lower bound
abline(v = x[max(which(xp < 0.1, arr.ind = T))])
##-- upper bound
abline(v = x[min(which(xp >= 0.9, arr.ind = T))])
text(x = lower - 0.6, y = 0.2, labels = paste0('lower\n', lower))
text(x = upper + 0.6, y = 0.2, labels = paste0('upper\n', upper))

k = 1
for(i in 1:nrow(chris)){
  text(x = chris$gene_cnlr_mean[i], y = k, 
       labels = chris$gene[i], 
       adj = 0.5, 
       cex = 0.7)
  k = k-0.06
}
title(main = paste0(sample, ' (10-90%)'))
dev.off()
write.table(x = chris, file = paste0('09_Cell_Line/', sample, '/', sample, '_IMPACT-like_CnLR.txt'), sep = '\t', row.names = F, quote = F)
rm(Mean, Sd, x, xp, snps, lower, upper, k, fourth, chris)


##-- FACETS
cval = 100
seed = 100
min_het = 15
snp_nbhd = 250

out = facetsSuite::run_facets(read_counts = countmatrix,
                              cval = cval,
                              dipLogR = NULL,
                              snp_nbhd = snp_nbhd,
                              seed = seed, 
                              genome = 'hg19', 
                              ndepth = 25)
out$dipLogR
i = cnlr_plot(facets_data = out, genome = 'hg19')
ii = valor_plot(facets_data = out, genome = 'hg19')
iii = icn_plot(facets_data = out, genome = 'hg19')
iv = cf_plot(facets_data = out, genome = 'hg19')

i / ii / iii / iv + plot_layout(heights = c(1,1,0.5,0.25))

qc = facets_fit_qc(facets_output = out)
qc

pass2 = i / ii / iii / iv + plot_layout(heights = c(1,1,0.5,0.25))
ggsave(filename = paste0('09_Cell_Line/', sample, '/', sample, '_facets.png'), plot = pass2, device = 'png', width = 12, height = 10)
saveRDS(object = out, file = '09_Cell_Line/FROZENPOOLEDNORMAL_IMPACT505_V2.rg.md.abra.printreads__s_TS543_14452/RDATA.rds')






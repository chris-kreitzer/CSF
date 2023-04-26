## Automated pipeline for CNA detection
## by using Facets other other technologies
## 
## start: 09/11/2022
## revision: 09/28/2022
## revision: 11/25/2022
## revision: 04/15/2023
## revision: 04/26/2023
## 
## chris-kreitzer

clean()
gc()
.rs.restartR()
setup(working.path = '~/Documents/MSKCC/11_CSF/')
library(patchwork)
library(dplyr)
library(data.table)
library(readr)
library(diptest)
library(mclust)
library(ggpubr)
source('~/Documents/GitHub/CSF/Scripts/UtilityFunctions.R')
source('~/Documents/GitHub/CSF/Scripts/FacetsPlot.R')
source('~/Documents/GitHub/CSF/Scripts/cf_Plot.R')
source('~/Documents/GitHub/CSF/Scripts/ReadSnpMatrix.R')
source('~/Documents/GitHub/CSF/Scripts/gene_closeup.R')
source('~/Documents/GitHub/CSF/Scripts/GMM.R')


database = readxl::read_excel('00_Data/Database_Chris_Sub_April12.xlsx')
csf = database[which(database$TYPE == 'CSF'), ]
csf = csf[which(csf$CSF_STATUS != 'Test Failure'), ]
files = list.files(path = '08_pileups/', full.names = T)
Nic_manifest = read.csv('06_Nic_Socci_r_001/Manifest_Nic.txt', sep = '\t')
GOIs = c('CDKN2A','CDK4','CDK6','PTEN','EGFR','PDGFRA',
         'KIT','KDR','MET','MDM2','MDM4','RB1','NF1',
         'TP53','FGF4', 'FGF19', 'PIK3CA', 'BRAF')


for(i in 1:nrow(csf)){
  if(any(grepl(pattern = csf$Sample.ID[i], files))){
    next
  } else {
    print(csf$Sample.ID[i])
  }
}





##----------------+
## First run
## countMatrix pre-check:
##----------------+
number = 160
sample = csf$Sample.ID[number]

countmatrix = readsnpmatrix(path = files[grep(pattern = sample, x = files)])
countmatrix = as.data.frame(countmatrix)

snps = facetsSuite::run_facets(read_counts = countmatrix,
                               cval = 100,
                               dipLogR = NULL,
                               snp_nbhd = 250,
                               seed = 100, 
                               genome = 'hg19', 
                               ndepth = 25)
het_snps = snps$snps[which(snps$snps$het == 1), ]
dim(het_snps)

dir.create(path = paste0('07_CSF_refit/', sample))
pdf(file = paste0('07_CSF_refit/', sample, '/', 'Normal_Het_Distribution.pdf'), width = 4.5, height = 4.5)
plot(density(het_snps$vafN[which(het_snps$het == 1)]), main = 'allele freq. het. SNPs\nmatched NORMAL')
dev.off()
norm_density = density(het_snps$vafN)
norm_density_max = norm_density$x[which.max(norm_density$y)]
ifelse(norm_density_max <= 0.4 | norm_density_max >= 0.6, print('STOP'), print('OK'))
pdf(file = paste0('07_CSF_refit/', sample, '/', 'TUMOR_Het_Distribution.pdf'), width = 4.5, height = 4.5)
plot(density(het_snps$vafT[which(het_snps$het == 1)]), main = 'allele freq. het. SNPs\nTUMOR')
dev.off()
rm(snps, het_snps)




##----------------+
## Second run:
## Ploidy, Purity, dipLogR determination
## (exclusively purity runs);
##----------------+
cval = 100
seed = 100
min_het = 15
snp_nbhd = 250

out = facetsSuite::run_facets(read_counts = countmatrix,
                              cval = cval,
                              dipLogR = 0.03,
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
system(command = paste0('mv ', files[grep(pattern = sample, x = files)], ' 07_CSF_refit/', sample, '/'))

pass2 = i / ii / iii / iv + plot_layout(heights = c(1,1,0.5,0.25))
ggsave(filename = paste0('07_CSF_refit/', sample, '/', sample, '_facets.png'), plot = pass2, device = 'png', width = 12, height = 10)

saveRDS(object = out, file = paste0('07_CSF_refit/', sample, '/', sample, '_second_pass.rds'))
out$dipLogR
rm(i, ii, iii, iv, pass2, out, cval, min_het)



##----------------+
## Third run:
## Gene-Level alterations if all the steps before
## are satisfied
## in-depth analyzes of CKDN2A, CDK6, and EGFR
##----------------+
cval = 50
seed = 100
min_het = 15
genome = 'hg19'
snp_nbhd = 100
diplogr = 0.03


fit = facetsSuite::run_facets(read_counts = countmatrix,
                              cval = 40,
                              dipLogR = diplogr,
                              snp_nbhd = snp_nbhd,
                              seed = seed, 
                              genome = 'hg19', 
                              ndepth = 25)

i = cnlr_plot(facets_data = fit, genome = 'hg19')
ii = valor_plot(facets_data = fit, genome = 'hg19')
iii = icn_plot(facets_data = fit, genome = 'hg19')
iv = cf_plot(facets_data = fit, genome = 'hg19')

i / ii / iii / iv + plot_layout(heights = c(1,1,0.5,0.25))


##-- gene level alteration
genes_all = facetsSuite::gene_level_changes(facets_output = fit, genome = 'hg19')
genes_all[which(genes_all$gene == 'EGFR'), c('chrom', 'median_cnlr_seg', 'tcn.em', 'lcn.em', 'cn_state', 'filter')]
genes_all[which(genes_all$gene == 'CDKN2A'), c('chrom', 'median_cnlr_seg', 'tcn.em', 'lcn.em', 'cn_state', 'filter')]
genes_all[which(genes_all$gene == 'CDK6'), c('chrom', 'median_cnlr_seg', 'tcn.em', 'lcn.em', 'cn_state', 'filter')]
genes_all = genes_all[which(genes_all$gene %in% GOIs), ]

write.table(x = genes_all, file = paste0('07_CSF_refit/', sample, '/', sample, '_gene_level_alteration.txt'), sep = '\t', row.names = F)
saveRDS(object = fit, file = paste0('07_CSF_refit/', sample, '/', sample, '_third_pass.rds'))

rm(fit, norm_density, i, ii, iii, iv)



##----------------+
## Fourth run: Engage with other
## methods: Nic Socci's calls and 
## the CnLR pipeline from my side
##----------------+
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

pdf(file = paste0('07_CSF_refit/', sample, '/', 'IMPACT-like.CnLR_distribution.pdf'), 
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
write.table(x = chris, file = paste0('07_CSF_refit/', sample, '/', sample, '_IMPACT-like_CnLR.txt'), sep = '\t', row.names = F, quote = F)

rm(Mean, Sd, x, xp, snps, lower, upper, k, fourth)



##----------------+
## check Nic Soccis calls
##----------------+

none




##----------------+
## Fifth run:
## Gaussian mixture model at selected genes
## with Cnlr only measures; CDKN2A, EGFR, CDK4, CDK6
##----------------+
dense_out = facetsSuite::run_facets(read_counts = countmatrix,
                                    cval = 100, 
                                    snp_nbhd = 50, 
                                    ndepth = 25, 
                                    seed = seed)
  

##-- investigate EGFR and CDKN2A manually
x = gene_closeup(data = dense_out, gene = 'CDKN2A')
x$plot
ggsave(filename = paste0('07_CSF_refit/', sample, '/', sample, '_CDKN2A_closeup.png'), plot = x$plot, 
       device = 'png', width = 10, height = 8)

sn = x$snps
hist(sn$cnlr, nclass = 100)
normality = ggqqplot(sn$cnlr, title = 'CnLR distribution chromosome 9p')
normality
ggsave(filename = paste0('07_CSF_refit/', sample, '/', sample, '_normality_9p.png'), plot = normality,
       device = 'png', width = 6, height = 6)


## MClust model; how many cluster fit the data
x.gmm = Mclust(sn$cnlr)
summary(x.gmm)
x.gmm$parameters$mean

vec1 = which(x.gmm$classification == 1, arr.ind = T)
jj = sn[vec1, ]
round(table(jj$gene)['color'][[1]] / sum(table(jj$gene))*100)

# none of the extreme cluster 1 refer to CDKN2A



##-- Gaussian mixture model; are there two components?
x = GMM(data = sn, components = 2)
x$plot
ggsave(filename = paste0('07_CSF_refit/', sample, '/', sample, '_CDKN2A_GMM.png'), plot = x$plot,
       device = 'png', width = 6, height = 6)

rm(x, jj, vec1, normality, sn, x.gmm)



##----------------+
## EGFR
##----------------+
x = gene_closeup(data = dense_out, gene = 'EGFR')
x$plot
ggsave(filename = paste0('07_CSF_refit/', sample, '/', sample, '_EGFR_closeup.png'), plot = x$plot, 
       device = 'png', width = 10, height = 8)

sn = x$snps

## MClust model; how many cluster fit the data
x.gmm = Mclust(sn$cnlr)
summary(x.gmm)
x.gmm$parameters$mean

vec1 = which(x.gmm$classification == 3, arr.ind = T)
jj = sn[vec1, ]
round(table(jj$gene)['color'][[1]] / sum(table(jj$gene))*100)

# none 

##-- Gaussian mixture model; are there two components?
x = GMM(data = sn, components = 3)
x$plot
ggsave(filename = paste0('07_CSF_refit/', sample, '/', sample, '_EGFR_GMM.png'), plot = x$plot,
       device = 'png', width = 6, height = 6)

rm(x, jj, vec1, sn, x.gmm)



sample_summary = data.frame(id = sample,
                            CNA_fit = 'pass',
                            Reason = '',
                            Notes = '',
                            Purity = qc$purity,
                            Ploidy = qc$ploidy,
                            FGA = qc$fga,
                            GMM = 'indication of CDKN2A and EGFR alteration',
                            Highlevel_CNA = c('EGFR_amplification', 'CDKN2A_deep_deletion'))

write.table(x = sample_summary, file = paste0('07_CSF_refit/', sample, '/', sample, '_summary.txt'), sep = '\t', row.names = F, quote = F)
rm(chris, countmatrix, dense_out, gene_out)










sample_summary = data.frame(id = sample,
                            CNA_fit = 'pass',
                            Reason = 'contamination; B-allele-frequency',
                            Notes = 'Careful. No purity,ploidy,fga,arm-level estimation possible',
                            Purity = NA,
                            Ploidy = NA,
                            GMM = 'no indication of CDKN2A alteration. EGFR homozygous deletion',
                            Highlevel_CNA = c('EGFR_amplification', 'CDKN2A_deep_deletion'))



## out
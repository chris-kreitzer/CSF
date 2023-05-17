setup(working.path = '~/Documents/MSKCC/12_CDKN2A/')
library(patchwork)
library(dplyr)
library(data.table)
library(readr)
library(mclust)
library(ggpubr)
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



sample = 'TCGA-41-3393-10A-01D-1353-08'
countmatrix = readsnpmatrix(path = 'TCGA-41-3393-10A-01D-1353-08_TCGA-41-3393-01A-01D-1353-08_C5o0bmTF_countsMerged___normal_tumor.dat.gz')
countmatrix = as.data.frame(countmatrix)



##-- Chris method
dir.create(path = paste0(sample))

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

pdf(file = paste0(sample, '/', 'IMPACT-like.CnLR_distribution.pdf'), 
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
write.table(x = chris, file = paste0(sample, '/', sample, '_IMPACT-like_CnLR.txt'), sep = '\t', row.names = F, quote = F)




##----------------+
## Facets run;
##----------------+
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
ggsave(filename = paste0(sample, '/', sample, '_cval100_facets.png'), plot = pass2, device = 'png', width = 12, height = 10)
saveRDS(object = out, file = paste0(sample, '/', sample, '_second_pass.rds'))

genes_all = facetsSuite::gene_level_changes(facets_output = out, genome = 'hg19')
genes_all[which(genes_all$gene == 'EGFR'), c('chrom', 'median_cnlr_seg', 'tcn.em', 'lcn.em', 'cn_state', 'filter')]
genes_all[which(genes_all$gene == 'CDKN2A'), c('chrom', 'median_cnlr_seg', 'tcn.em', 'lcn.em', 'cn_state', 'filter')]



##----------------+
## GMM based
##----------------+
dense_out = fourth

##----------------+
## CDKN2A
##----------------+
x = gene_closeup(data = dense_out, gene = 'CDKN2A')
arm = x$plot
exon = exon_mapping(data = x$snps)

cdkn2a = arm / exon
diplogr = dense_out$dipLogR
cdkn2a
ggsave(filename = paste0(sample, '/', sample, '_CDKN2A_closeup.png'), plot = cdkn2a, 
       device = 'png', width = 12, height = 10)





sn = x$snps
## MClust model; how many cluster fit the data
x.gmm = Mclust(sn$cnlr)
summary(x.gmm)
x.gmm$parameters$mean


total = table(sn$gene)[['color']][1]
gmm_out_all = data.frame()
for(i in 1:length(x.gmm$parameters$mean)){
  cl_m = round(x.gmm$parameters$mean[i][[1]], 2)
  vec = which(x.gmm$classification == i, arr.ind = T)
  jj = sn[vec, ]
  cluster = round(table(jj$gene)['color'][[1]] / total, 2) * 100
  gmm_out = data.frame(Cluster = i,
                       cl_m = cl_m,
                       cluster_proportion = cluster,
                       gene = 'CDKN2A')
  gmm_out_all = rbind(gmm_out_all, gmm_out)
  
}
gmm_out_all


##-- Gaussian mixture model; are there two components?
x = GMM(data = sn, components = length(x.gmm$parameters$mean))
plot = x$plot
diplogr = round(diplogr, 3)
plot = plot + annotate(geom = 'text',
                       x = -1.5,
                       y = 0.5, 
                       label = paste0('DipLogR of sample: ', diplogr, '\n',
                                      'Total probes: ', total, '\n',
                                      'cluster 1: ', gmm_out_all$cl_m[which(gmm_out_all$Cluster == 1)], '; ', gmm_out_all$cluster_proportion[which(gmm_out_all$Cluster == 1)], '% of probes\n',
                                      'cluster 2: ', gmm_out_all$cl_m[which(gmm_out_all$Cluster == 2)], '; ', gmm_out_all$cluster_proportion[which(gmm_out_all$Cluster == 2)], '% of probes\n',
                                      'cluster 3: ', gmm_out_all$cl_m[which(gmm_out_all$Cluster == 3)], '; ', gmm_out_all$cluster_proportion[which(gmm_out_all$Cluster == 3)], '% of probes\n')) +
  labs(title = 'CDKN2A; Fraction of probes falling within cluster')

plot
ggsave(filename = paste0(sample, '/', sample, '_CDKN2A_GMM.png'), plot = plot,
       device = 'png', width = 9, height = 6)

genes_all[which(genes_all$gene == 'CDKN2A'), c('chrom', 'median_cnlr_seg', 'tcn.em', 'lcn.em', 'cn_state', 'filter')]

rm(x, jj, sn, x.gmm, plot, gmm_out, gmm_out_all)






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
source('~/Documents/GitHub/CSF/Scripts/ExonicMapping_CDKN2A.R')


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
## Normal/Tumor pairs
##----------------+
number = 166
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
## Second run
## check Nic Soccis calls
##----------------+



#----------------+
## Third run: 
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
## Fourth run:
## FACETS: Ploidy, Purity, dipLogR determination
## (exclusively purity runs);
## broad run
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
system(command = paste0('mv ', files[grep(pattern = sample, x = files)], ' 07_CSF_refit/', sample, '/'))

pass2 = i / ii / iii / iv + plot_layout(heights = c(1,1,0.5,0.25))
ggsave(filename = paste0('07_CSF_refit/', sample, '/', sample, '_facets.png'), plot = pass2, device = 'png', width = 12, height = 10)

saveRDS(object = out, file = paste0('07_CSF_refit/', sample, '/', sample, '_second_pass.rds'))
out$dipLogR
rm(i, ii, iii, iv, pass2, out, cval, min_het)
genes_broad = facetsSuite::gene_level_changes(facets_output = out, genome = 'hg19')



##----------------+
## Fifth run:
## FACETS- hisens run
## Gene-Level alterations if all the steps before
## are satisfied
## in-depth analyzes of CKDN2A, CDK6, and EGFR
##----------------+
cval = 50
seed = 100
min_het = 15
genome = 'hg19'
snp_nbhd = 100
diplogr = -0.07208348


fit = facetsSuite::run_facets(read_counts = countmatrix,
                              cval = 50,
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
genes_all[which(genes_all$gene == 'CDK4'), c('chrom', 'median_cnlr_seg', 'tcn.em', 'lcn.em', 'cn_state', 'filter')]
genes_all[which(genes_all$gene == 'PTEN'), c('chrom', 'median_cnlr_seg', 'tcn.em', 'lcn.em', 'cn_state', 'filter')]
genes_all = genes_all[which(genes_all$gene %in% GOIs), ]

#write.table(x = genes_all, file = paste0('07_CSF_refit/', sample, '/', sample, '_gene_level_alteration.txt'), sep = '\t', row.names = F)
#saveRDS(object = fit, file = paste0('07_CSF_refit/', sample, '/', sample, '_third_pass.rds'))

#rm(fit, norm_density, i, ii, iii, iv)
genes_all[,c('gene', 'chrom', 'tcn.em', 'cn_state', 'filter')]


##--- choose which FACETS call to use (broad / hisens run)



##----------------+
## Sixth run:
## Gaussian mixture model at selected genes
## with Cnlr only measures; CDKN2A, EGFR, CDK4, CDK6
##----------------+
countmatrix = readsnpmatrix(path = '01_countmatrices/countsMerged____P-0005701-T01-IM5_P-0005701-N01-IM5.dat.gz/countsMerged____P-0005701-T01-IM5_P-0005701-N01-IM5.dat.gz')
dense_out = facetsSuite::run_facets(read_counts = countmatrix,
                                    cval = 100, 
                                    snp_nbhd = 50, 
                                    ndepth = 25, 
                                    seed = seed)
  

##----------------+
## CDKN2A
##----------------+
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



##----------------+
## exonic structure CDKN2A
##----------------+
cd = sn[which(sn$maploc >= 21967752 & sn$maploc <= 21995300), ]
cd$seq = 1:nrow(cd)
cd$seq = cd$seq * -1


exonic = exonic_structure(gene = 'CDKN2A', type = 'exon')
exonic = exonic[which(exonic$level  == 1), ]

View(exonic)

exon3: 21967751-21968241
exon2: 21970901-21971207
exon1a

ggplot(cd) +
  geom_rect(aes(xmin = 21970000, xmax = 21975000, ymin = 0.5, ymax = 1.5)) +
  geom_point(aes(x = maploc, y = seq, color = cnlr), size = 2) +
  scale_colour_gradient2(
    low = "darkblue",
    mid = "white",
    high = "red",
    midpoint = -1.7,
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour"
  ) +

  theme_bw() +
  theme(aspect.ratio = 0.5)
  




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
plot = plot + annotate(geom = 'text',
                x = -2,
                y = 0.5, 
                label = paste0('DipLogR of sample: ', diplogr, '\n',
                               'Total probes: ', total, '\n',
                               'cluster 1: ', gmm_out_all$cl_m[which(gmm_out_all$Cluster == 1)], '; ', gmm_out_all$cluster_proportion[which(gmm_out_all$Cluster == 1)], '% of probes\n',
                               'cluster 2: ', gmm_out_all$cl_m[which(gmm_out_all$Cluster == 2)], '; ', gmm_out_all$cluster_proportion[which(gmm_out_all$Cluster == 2)], '% of probes\n',
                               'cluster 3: ', gmm_out_all$cl_m[which(gmm_out_all$Cluster == 3)], '; ', gmm_out_all$cluster_proportion[which(gmm_out_all$Cluster == 3)], '% of probes\n')) +
  labs(title = 'CDKN2A; Fraction of probes falling within cluster')

plot
ggsave(filename = paste0('07_CSF_refit/', sample, '/', sample, '_CDKN2A_GMM.png'), plot = plot,
       device = 'png', width = 9, height = 6)

genes_all[which(genes_all$gene == 'CDKN2A'), c('chrom', 'median_cnlr_seg', 'tcn.em', 'lcn.em', 'cn_state', 'filter')]
genes_broad[which(genes_broad$gene == 'CDKN2A'), c('chrom', 'median_cnlr_seg', 'tcn.em', 'lcn.em', 'cn_state', 'filter')]

rm(x, jj, vec1, normality, sn, x.gmm, plot, gmm_out, gmm_out_all)



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
                       gene = 'EGFR')
  gmm_out_all = rbind(gmm_out_all, gmm_out)
  
}
gmm_out_all



##-- Gaussian mixture model; are there two components?
x = GMM(data = sn, components = length(x.gmm$parameters$mean))
plot = x$plot
plot = plot + annotate(geom = 'text',
                       x = -2,
                       y = 0.3, 
                       label = paste0('DipLogR of sample: ', diplogr, '\n',
                                      'Total probes: ', total, '\n',
                                      'cluster 1: ', gmm_out_all$cl_m[which(gmm_out_all$Cluster == 1)], '; ', gmm_out_all$cluster_proportion[which(gmm_out_all$Cluster == 1)], '% of probes\n',
                                      'cluster 2: ', gmm_out_all$cl_m[which(gmm_out_all$Cluster == 2)], '; ', gmm_out_all$cluster_proportion[which(gmm_out_all$Cluster == 2)], '% of probes\n',
                                      'cluster 3: ', gmm_out_all$cl_m[which(gmm_out_all$Cluster == 3)], '; ', gmm_out_all$cluster_proportion[which(gmm_out_all$Cluster == 3)], '% of probes\n')) +
  labs(title = 'EGFR; Fraction of probes falling within cluster')

plot
ggsave(filename = paste0('07_CSF_refit/', sample, '/', sample, '_EGFR_GMM.png'), plot = plot,
       device = 'png', width = 9, height = 6)

genes_all[which(genes_all$gene == 'EGFR'), c('chrom', 'median_cnlr_seg', 'tcn.em', 'lcn.em', 'cn_state', 'filter')]
genes_broad[which(genes_broad$gene == 'EGFR'), c('chrom', 'median_cnlr_seg', 'tcn.em', 'lcn.em', 'cn_state', 'filter')]

rm(x, jj, vec1, vec2, vec3, cl2, cl1, cl3, jj, total, sn, x.gmm, plot, gmm_out, gmm_out_all)



##----------------+
## CDK4
##----------------+
x = gene_closeup(data = dense_out, gene = 'CDK4')
x$plot
ggsave(filename = paste0('07_CSF_refit/', sample, '/', sample, '_CDK4_closeup.png'), plot = x$plot, 
       device = 'png', width = 10, height = 8)

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
                       gene = 'CDK4')
  gmm_out_all = rbind(gmm_out_all, gmm_out)
  
}
gmm_out_all



##-- Gaussian mixture model; are there two components?
x = GMM(data = sn, components = length(x.gmm$parameters$mean))
plot = x$plot
plot = plot + annotate(geom = 'text',
                       x = -2,
                       y = 0.3, 
                       label = paste0('DipLogR of sample: ', diplogr, '\n',
                                      'Total probes: ', total, '\n',
                                      'cluster 1: ', gmm_out_all$cl_m[which(gmm_out_all$Cluster == 1)], '; ', gmm_out_all$cluster_proportion[which(gmm_out_all$Cluster == 1)], '% of probes\n',
                                      'cluster 2: ', gmm_out_all$cl_m[which(gmm_out_all$Cluster == 2)], '; ', gmm_out_all$cluster_proportion[which(gmm_out_all$Cluster == 2)], '% of probes\n',
                                      'cluster 3: ', gmm_out_all$cl_m[which(gmm_out_all$Cluster == 3)], '; ', gmm_out_all$cluster_proportion[which(gmm_out_all$Cluster == 3)], '% of probes\n')) +
  labs(title = 'CDK4; Fraction of probes falling within cluster')

plot
ggsave(filename = paste0('07_CSF_refit/', sample, '/', sample, '_CDK4_GMM.png'), plot = plot,
       device = 'png', width = 9, height = 6)
genes_all[which(genes_all$gene == 'CDK4'), c('chrom', 'median_cnlr_seg', 'tcn.em', 'lcn.em', 'cn_state', 'filter')]
genes_broad[which(genes_broad$gene == 'CDK4'), c('chrom', 'median_cnlr_seg', 'tcn.em', 'lcn.em', 'cn_state', 'filter')]

rm(x, vec1, vec2, vec3, cl2, cl1, cl3, jj, total, sn, x.gmm, plot)



##----------------+
## CDK6
##----------------+
x = gene_closeup(data = dense_out, gene = 'CDK6')
x$plot
ggsave(filename = paste0('07_CSF_refit/', sample, '/', sample, '_CDK6_closeup.png'), plot = x$plot, 
       device = 'png', width = 10, height = 8)

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
                       gene = 'CDK6')
  gmm_out_all = rbind(gmm_out_all, gmm_out)
  
}
gmm_out_all



##-- Gaussian mixture model; are there two components?
x = GMM(data = sn, components = length(x.gmm$parameters$mean))
plot = x$plot
plot = plot + annotate(geom = 'text',
                       x = -2,
                       y = 0.3, 
                       label = paste0('DipLogR of sample: ', diplogr, '\n',
                                      'Total probes: ', total, '\n',
                                      'cluster 1: ', gmm_out_all$cl_m[which(gmm_out_all$Cluster == 1)], '; ', gmm_out_all$cluster_proportion[which(gmm_out_all$Cluster == 1)], '% of probes\n',
                                      'cluster 2: ', gmm_out_all$cl_m[which(gmm_out_all$Cluster == 2)], '; ', gmm_out_all$cluster_proportion[which(gmm_out_all$Cluster == 2)], '% of probes\n',
                                      'cluster 3: ', gmm_out_all$cl_m[which(gmm_out_all$Cluster == 3)], '; ', gmm_out_all$cluster_proportion[which(gmm_out_all$Cluster == 3)], '% of probes\n')) +
  labs(title = 'CDK6; Fraction of probes falling within cluster')

plot
ggsave(filename = paste0('07_CSF_refit/', sample, '/', sample, '_CDK6_GMM.png'), plot = plot,
       device = 'png', width = 9, height = 6)
genes_all[which(genes_all$gene == 'CDK6'), c('chrom', 'median_cnlr_seg', 'tcn.em', 'lcn.em', 'cn_state', 'filter')]
genes_broad[which(genes_broad$gene == 'CDK6'), c('chrom', 'median_cnlr_seg', 'tcn.em', 'lcn.em', 'cn_state', 'filter')]


rm(x, vec1, vec2, vec3, cl2, cl1, cl3, jj, total, sn, x.gmm, plot)




##----------------+
## PTEN
##----------------+
x = gene_closeup(data = dense_out, gene = 'PTEN')
x$plot
ggsave(filename = paste0('07_CSF_refit/', sample, '/', sample, '_PTEN_closeup.png'), plot = x$plot, 
       device = 'png', width = 10, height = 8)

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
                       gene = 'PTEN')
  gmm_out_all = rbind(gmm_out_all, gmm_out)
  
}
gmm_out_all



##-- Gaussian mixture model; are there two components?
x = GMM(data = sn, components = length(x.gmm$parameters$mean))
plot = x$plot
plot = plot + annotate(geom = 'text',
                       x = -2,
                       y = 0.3, 
                       label = paste0('DipLogR of sample: ', diplogr, '\n',
                                      'Total probes: ', total, '\n',
                                      'cluster 1: ', gmm_out_all$cl_m[which(gmm_out_all$Cluster == 1)], '; ', gmm_out_all$cluster_proportion[which(gmm_out_all$Cluster == 1)], '% of probes\n',
                                      'cluster 2: ', gmm_out_all$cl_m[which(gmm_out_all$Cluster == 2)], '; ', gmm_out_all$cluster_proportion[which(gmm_out_all$Cluster == 2)], '% of probes\n',
                                      'cluster 3: ', gmm_out_all$cl_m[which(gmm_out_all$Cluster == 3)], '; ', gmm_out_all$cluster_proportion[which(gmm_out_all$Cluster == 3)], '% of probes\n')) +
  labs(title = 'PTEN; Fraction of probes falling within cluster')

plot
ggsave(filename = paste0('07_CSF_refit/', sample, '/', sample, '_PTEN_GMM.png'), plot = plot,
       device = 'png', width = 9, height = 6)
genes_all[which(genes_all$gene == 'PTEN'), c('chrom', 'median_cnlr_seg', 'tcn.em', 'lcn.em', 'cn_state', 'filter')]
genes_broad[which(genes_broad$gene == 'PTEN'), c('chrom', 'median_cnlr_seg', 'tcn.em', 'lcn.em', 'cn_state', 'filter')]


rm(x, vec1, vec2, vec3, cl2, cl1, cl3, jj, total, sn, x.gmm, plot)



## continue herre

##----------------+
## MDM2
##----------------+
x = gene_closeup(data = dense_out, gene = 'MDM2')
x$plot
ggsave(filename = paste0('07_CSF_refit/', sample, '/', sample, '_MDM2_closeup.png'), plot = x$plot, 
       device = 'png', width = 10, height = 8)

sn = x$snps

## MClust model; how many cluster fit the data
## MClust model; how many cluster fit the data
x.gmm = Mclust(sn$cnlr)
summary(x.gmm)
x.gmm$parameters$mean

cl1_m = round(x.gmm$parameters$mean[1][[1]], 2)
cl2_m = round(x.gmm$parameters$mean[2][[1]], 2)
cl3_m = round(x.gmm$parameters$mean[3][[1]], 2)

diplogr
total = table(sn$gene)[['color']][1]
total

vec1 = which(x.gmm$classification == 1, arr.ind = T)
jj = sn[vec1, ]
cl1 = round(table(jj$gene)['color'][[1]] / total, 2) * 100
cl1

vec2 = which(x.gmm$classification == 2, arr.ind = T)
jj = sn[vec2, ]
cl2 = round(table(jj$gene)['color'][[1]] / total, 2) * 100
cl2 

vec3 = which(x.gmm$classification == 4, arr.ind = T)
jj = sn[vec3, ]
cl3 = round(table(jj$gene)['color'][[1]] / total, 2) * 100
cl3


##-- Gaussian mixture model; are there two components?
x = GMM(data = sn, components = 2)
plot = x$plot
plot = plot + annotate(geom = 'text',
                       x = -2,
                       y = 0.3, 
                       label = paste0('DipLogR of sample: ', diplogr, '\n',
                                      'cluster 1: ', cl1_m, '; ', cl1, '% of probes\n',
                                      'cluster 2: ', cl2_m, '; ', cl2, '% of probes\n',
                                      'cluster 3: ', cl3_m, '; ', cl3, '% of probes\n')) +
  labs(title = 'MDM2; Fraction of probes falling within cluster')

plot
ggsave(filename = paste0('07_CSF_refit/', sample, '/', sample, '_MDM2_GMM.png'), plot = plot,
       device = 'png', width = 9, height = 6)
genes_all[which(genes_all$gene == 'MDM2'), c('chrom', 'median_cnlr_seg', 'tcn.em', 'lcn.em', 'cn_state', 'filter')]
genes_broad[which(genes_broad$gene == 'MDM2'), c('chrom', 'median_cnlr_seg', 'tcn.em', 'lcn.em', 'cn_state', 'filter')]


rm(x, vec1, vec2, vec3, cl2, cl1, cl3, jj, total, sn, x.gmm, plot)



##----------------+
## KIT
##----------------+
x = gene_closeup(data = dense_out, gene = 'KIT')
x$plot
ggsave(filename = paste0('07_CSF_refit/', sample, '/', sample, '_KIT_closeup.png'), plot = x$plot, 
       device = 'png', width = 10, height = 8)

sn = x$snps

## MClust model; how many cluster fit the data
## MClust model; how many cluster fit the data
x.gmm = Mclust(sn$cnlr)
summary(x.gmm)
x.gmm$parameters$mean

cl1_m = round(x.gmm$parameters$mean[1][[1]], 2)
cl2_m = round(x.gmm$parameters$mean[2][[1]], 2)
cl3_m = round(x.gmm$parameters$mean[3][[1]], 2)

diplogr
total = table(sn$gene)[['color']][1]
total

vec1 = which(x.gmm$classification == 1, arr.ind = T)
jj = sn[vec1, ]
cl1 = round(table(jj$gene)['color'][[1]] / total, 2) * 100
cl1

vec2 = which(x.gmm$classification == 2, arr.ind = T)
jj = sn[vec2, ]
cl2 = round(table(jj$gene)['color'][[1]] / total, 2) * 100
cl2 

vec3 = which(x.gmm$classification == 3, arr.ind = T)
jj = sn[vec3, ]
cl3 = round(table(jj$gene)['color'][[1]] / total, 2) * 100
cl3

vec4 = which(x.gmm$classification == 4, arr.ind = T)
jj = sn[vec4, ]
cl4 = round(table(jj$gene)['color'][[1]] / total, 2) * 100
cl4

##-- Gaussian mixture model; are there two components?
x = GMM(data = sn, components = 3)
plot = x$plot
plot = plot + annotate(geom = 'text',
                       x = -2,
                       y = 0.3, 
                       label = paste0('DipLogR of sample: ', diplogr, '\n',
                                      'cluster 1: ', cl1_m, '; ', cl1, '% of probes\n',
                                      'cluster 2: ', cl2_m, '; ', cl2, '% of probes\n',
                                      'cluster 3: ', cl3_m, '; ', cl3, '% of probes\n')) +
  labs(title = 'KIT; Fraction of probes falling within cluster')

plot
ggsave(filename = paste0('07_CSF_refit/', sample, '/', sample, '_KIT_GMM.png'), plot = plot,
       device = 'png', width = 9, height = 6)

rm(x, vec1, vec2, vec3, cl2, cl1, cl3, jj, total, sn, x.gmm, plot)



##----------------+
## MET
##----------------+
x = gene_closeup(data = dense_out, gene = 'MET')
x$plot
ggsave(filename = paste0('07_CSF_refit/', sample, '/', sample, '_MET_closeup.png'), plot = x$plot, 
       device = 'png', width = 10, height = 8)

sn = x$snps

## MClust model; how many cluster fit the data
## MClust model; how many cluster fit the data
x.gmm = Mclust(sn$cnlr)
summary(x.gmm)
x.gmm$parameters$mean

cl1_m = round(x.gmm$parameters$mean[1][[1]], 2)
cl2_m = round(x.gmm$parameters$mean[2][[1]], 2)
cl3_m = round(x.gmm$parameters$mean[3][[1]], 2)

diplogr
total = table(sn$gene)[['color']][1]
total

vec1 = which(x.gmm$classification == 1, arr.ind = T)
jj = sn[vec1, ]
cl1 = round(table(jj$gene)['color'][[1]] / total, 2) * 100
cl1

vec2 = which(x.gmm$classification == 2, arr.ind = T)
jj = sn[vec2, ]
cl2 = round(table(jj$gene)['color'][[1]] / total, 2) * 100
cl2 

vec3 = which(x.gmm$classification == 4, arr.ind = T)
jj = sn[vec3, ]
cl3 = round(table(jj$gene)['color'][[1]] / total, 2) * 100
cl3


##-- Gaussian mixture model; are there two components?
x = GMM(data = sn, components = 3)
plot = x$plot
plot = plot + annotate(geom = 'text',
                       x = -2,
                       y = 0.3, 
                       label = paste0('DipLogR of sample: ', diplogr, '\n',
                                      'cluster 1: ', cl1_m, '; ', cl1, '% of probes\n',
                                      'cluster 2: ', cl2_m, '; ', cl2, '% of probes\n',
                                      'cluster 3: ', cl3_m, '; ', cl3, '% of probes\n')) +
  labs(title = 'MET; Fraction of probes falling within cluster')

plot
ggsave(filename = paste0('07_CSF_refit/', sample, '/', sample, '_MET_GMM.png'), plot = plot,
       device = 'png', width = 9, height = 6)

rm(x, vec1, vec2, vec3, cl2, cl1, cl3, jj, total, sn, x.gmm, plot)



##----------------+
## PDGFRA
##----------------+
x = gene_closeup(data = dense_out, gene = 'PDGFRA')
x$plot
ggsave(filename = paste0('07_CSF_refit/', sample, '/', sample, '_PDGFRA_closeup.png'), plot = x$plot, 
       device = 'png', width = 10, height = 8)

sn = x$snps

## MClust model; how many cluster fit the data
## MClust model; how many cluster fit the data
x.gmm = Mclust(sn$cnlr)
summary(x.gmm)
x.gmm$parameters$mean

cl1_m = round(x.gmm$parameters$mean[1][[1]], 2)
cl2_m = round(x.gmm$parameters$mean[2][[1]], 2)
cl3_m = round(x.gmm$parameters$mean[3][[1]], 2)

diplogr
total = table(sn$gene)[['color']][1]
total

vec1 = which(x.gmm$classification == 1, arr.ind = T)
jj = sn[vec1, ]
cl1 = round(table(jj$gene)['color'][[1]] / total, 2) * 100
cl1

vec2 = which(x.gmm$classification == 2, arr.ind = T)
jj = sn[vec2, ]
cl2 = round(table(jj$gene)['color'][[1]] / total, 2) * 100
cl2 

vec3 = which(x.gmm$classification == 4, arr.ind = T)
jj = sn[vec3, ]
cl3 = round(table(jj$gene)['color'][[1]] / total, 2) * 100
cl3


##-- Gaussian mixture model; are there two components?
x = GMM(data = sn, components = 3)
plot = x$plot
plot = plot + annotate(geom = 'text',
                       x = -2,
                       y = 0.3, 
                       label = paste0('DipLogR of sample: ', diplogr, '\n',
                                      'cluster 1: ', cl1_m, '; ', cl1, '% of probes\n',
                                      'cluster 2: ', cl2_m, '; ', cl2, '% of probes\n',
                                      'cluster 3: ', cl3_m, '; ', cl3, '% of probes\n')) +
  labs(title = 'PDGFRA; Fraction of probes falling within cluster')

plot
ggsave(filename = paste0('07_CSF_refit/', sample, '/', sample, '_PDGFRA_GMM.png'), plot = plot,
       device = 'png', width = 9, height = 6)

rm(x, vec1, vec2, vec3, cl2, cl1, cl3, jj, total, sn, x.gmm, plot)



##----------------+
## RB1
##----------------+
x = gene_closeup(data = dense_out, gene = 'RB1')
x$plot
ggsave(filename = paste0('07_CSF_refit/', sample, '/', sample, '_RB1_closeup.png'), plot = x$plot, 
       device = 'png', width = 10, height = 8)

sn = x$snps

## MClust model; how many cluster fit the data
## MClust model; how many cluster fit the data
x.gmm = Mclust(sn$cnlr)
summary(x.gmm)
x.gmm$parameters$mean

cl1_m = round(x.gmm$parameters$mean[1][[1]], 2)
cl2_m = round(x.gmm$parameters$mean[2][[1]], 2)
cl3_m = round(x.gmm$parameters$mean[3][[1]], 2)

diplogr
total = table(sn$gene)[['color']][1]
total

vec1 = which(x.gmm$classification == 1, arr.ind = T)
jj = sn[vec1, ]
cl1 = round(table(jj$gene)['color'][[1]] / total, 2) * 100
cl1

vec2 = which(x.gmm$classification == 2, arr.ind = T)
jj = sn[vec2, ]
cl2 = round(table(jj$gene)['color'][[1]] / total, 2) * 100
cl2 

vec3 = which(x.gmm$classification == 4, arr.ind = T)
jj = sn[vec3, ]
cl3 = round(table(jj$gene)['color'][[1]] / total, 2) * 100
cl3


##-- Gaussian mixture model; are there two components?
x = GMM(data = sn, components = 3)
plot = x$plot
plot = plot + annotate(geom = 'text',
                       x = -2,
                       y = 0.3, 
                       label = paste0('DipLogR of sample: ', diplogr, '\n',
                                      'cluster 1: ', cl1_m, '; ', cl1, '% of probes\n',
                                      'cluster 2: ', cl2_m, '; ', cl2, '% of probes\n',
                                      'cluster 3: ', cl3_m, '; ', cl3, '% of probes\n')) +
  labs(title = 'RB1; Fraction of probes falling within cluster')

plot
ggsave(filename = paste0('07_CSF_refit/', sample, '/', sample, '_RB1_GMM.png'), plot = plot,
       device = 'png', width = 9, height = 6)

rm(x, vec1, vec2, vec3, cl2, cl1, cl3, jj, total, sn, x.gmm, plot)



##----------------+
## KDR
##----------------+
x = gene_closeup(data = dense_out, gene = 'KDR')
x$plot
ggsave(filename = paste0('07_CSF_refit/', sample, '/', sample, 'KDR_closeup.png'), plot = x$plot, 
       device = 'png', width = 10, height = 8)

sn = x$snps

## MClust model; how many cluster fit the data
## MClust model; how many cluster fit the data
x.gmm = Mclust(sn$cnlr)
summary(x.gmm)
x.gmm$parameters$mean

cl1_m = round(x.gmm$parameters$mean[1][[1]], 2)
cl2_m = round(x.gmm$parameters$mean[2][[1]], 2)
cl3_m = round(x.gmm$parameters$mean[3][[1]], 2)

diplogr
total = table(sn$gene)[['color']][1]
total

vec1 = which(x.gmm$classification == 1, arr.ind = T)
jj = sn[vec1, ]
cl1 = round(table(jj$gene)['color'][[1]] / total, 2) * 100
cl1

vec2 = which(x.gmm$classification == 2, arr.ind = T)
jj = sn[vec2, ]
cl2 = round(table(jj$gene)['color'][[1]] / total, 2) * 100
cl2 

vec3 = which(x.gmm$classification == 4, arr.ind = T)
jj = sn[vec3, ]
cl3 = round(table(jj$gene)['color'][[1]] / total, 2) * 100
cl3


##-- Gaussian mixture model; are there two components?
x = GMM(data = sn, components = 3)
plot = x$plot
plot = plot + annotate(geom = 'text',
                       x = -2,
                       y = 0.3, 
                       label = paste0('DipLogR of sample: ', diplogr, '\n',
                                      'cluster 1: ', cl1_m, '; ', cl1, '% of probes\n',
                                      'cluster 2: ', cl2_m, '; ', cl2, '% of probes\n',
                                      'cluster 3: ', cl3_m, '; ', cl3, '% of probes\n')) +
  labs(title = 'KDR; Fraction of probes falling within cluster')

plot
ggsave(filename = paste0('07_CSF_refit/', sample, '/', sample, '_KDR_GMM.png'), plot = plot,
       device = 'png', width = 9, height = 6)

rm(x, vec1, vec2, vec3, cl2, cl1, cl3, jj, total, sn, x.gmm, plot)


##----------------+
## summarize the output
##----------------+
sample_summary = data.frame(id = sample,
                            CNA_fit = 'Manual',
                            Notes = 'GMM based. No purity, ploidy.',
                            Clonality_analysis = 'No',
                            Purity = NA,
                            Ploidy = NA,
                            FGA = NA,
                            CDKN2A = 'Diploid',
                            CDK4 = 'Diploid',
                            EGFR = 'Loss',
                            CDK6 = 'Loss',
                            PTEN = 'Diploid',
                            KIT = 'Loss',
                            MDM2 = 'Diploid',
                            MET = 'Loss',
                            RB1 = 'Diploid',
                            KDR = 'Loss',
                            Highlevel_CNA = NA)


write.table(x = sample_summary, file = paste0('07_CSF_refit/', sample, '/', sample, '_summary.txt'), sep = '\t', row.names = F, quote = F)
rm(chris, countmatrix, dense_out, gene_out, qc, sample_summary, position)
dev.off()


## out

setup(working.path = '~/Documents/MSKCC/11_CSF/')
library(patchwork)
library(dplyr)
library(data.table)
library(readr)
library(ggpubr)
library(cowplot)

source('~/Documents/GitHub/CSF/Scripts/UtilityFunctions.R')
source('~/Documents/GitHub/CSF/Scripts/FacetsPlot.R')
source('~/Documents/GitHub/CSF/Scripts/cf_Plot.R')
source('~/Documents/GitHub/CSF/Scripts/ReadSnpMatrix.R')

countmatrix = readsnpmatrix('01_countmatrices/countsMerged____P-0005701-T01-IM5_P-0005701-N01-IM5.dat.gz/countsMerged____P-0005701-T01-IM5_P-0005701-N01-IM5.dat.gz')
countmatrix = as.data.frame(countmatrix)



out = facetsSuite::run_facets(read_counts = countmatrix,
                              snp_nbhd = 50)


i = cnlr_plot(facets_data = out, genome = 'hg19')
ii = valor_plot(facets_data = out, genome = 'hg19')
iii = icn_plot(facets_data = out, genome = 'hg19')
iv = cf_plot(facets_data = out, genome = 'hg19')

i / ii / iii / iv + plot_layout(heights = c(1,1,0.5,0.25))


snps = out$snps
snps = snps[which(snps$chrom == 9), ]
hist(snps$cnlr, nclass = 100, xlab = 'CnLR', main = 'Chromosome 9 of\nP-0005701-T01-IM5')
ggqqplot(snps$cnlr)

##' 9p arm

p_snps = snps[which(snps$maploc < 48000000), ]
q_snps = snps[which(snps$maploc > 48000000), ]


hist(p_snps$cnlr, nclass = 100, 
     xlab = 'CnLR', 
     main = 'Chromosome-arm 9p of\nP-0005701-T01-IM5')
hist(q_snps$cnlr, nclass = 100, 
     xlab = 'CnLR', 
     main = 'Chromosome-arm 9q of\nP-0005701-T01-IM5')


##' visual exploration 
p_snps$order = seq(1, nrow(p_snps), 1)
p_snps$cdkn2a = ifelse(p_snps$maploc >= 21967751 & p_snps$maploc <= 21995300, 'color', 'no_color')

ggplot(p_snps, aes(x = order, y = cnlr, color = cdkn2a)) +
  geom_point(position = position_dodge(width = 0.1), size = 1) +
  scale_color_manual(values = c('color' = '#487EFB', 
                               'no_color' = 'grey35'),
                     name = '', labels = c('CDKN2A', 'other')) +
  scale_y_continuous(expand = c(0.01, 1), limits = c(-3, 2)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  panel_border(size = 2, color = 'black') +
  labs(x = 'Chromosome 9p', y = 'CnLR')



##' Gaussian mixture model - indication of subpopulation within broader loss
library(nortest)
ad.test(x = p_snps$cnlr)
cvm.test(x = snps$cnlr)
library(diptest)
dip.test(x = snps$cnlr, simulate.p.value = FALSE, B = 2000)

nn = rnorm(n = 100, mean = 2, sd = 1)
uu = rnorm(n = 100, mean = -2, sd = 1)
xx = c(nn, uu)

dip.test(x = xx)
hist(xx, nclass = 100)

library(mclust)
x.gmm = Mclust(p_snps$cnlr)
summary(x.gmm)
x.gmm$parameters
hist(p_snps$cnlr, nclass = 150)
abline(v = x.gmm$parameters$mean[1])
abline(v = x.gmm$parameters$mean[2])



mixmdl = mixtools::normalmixEM(x = p_snps$cnlr, k = 3, mu = c(-1,0, 1), ECM = F)
summary(mixmdl)














km.fromscratch <- function(X, k){
  p <- ncol(X)  # number of parameters
  n <- nrow(X)  # number of observations
  Delta <- 1; iter <- 0; itermax <- 30
  while(Delta > 1e-4 && iter <= itermax){
    # initiation
    if(iter == 0){
      centroid <- X[sample(nrow(X), k),]
      centroid_mem <- centroid
    }
    
    # equivalent to E-step
    d <- sapply(1:k, function(c) sapply(1:n, 
                                        function(i) sum((centroid[c,] - X[i,])^2) ))
    cluster <- apply(d, 1, which.min)
    
    # equivalent to M-step
    centroid <- t(sapply(1:k, function(c) 
      apply(X[cluster == c,], 2, mean)))
    
    Delta <- sum((centroid - centroid_mem)^2)
    iter <- iter + 1; centroid_mem <- centroid
  }
  return(list(centroid = centroid, cluster = cluster))
}
# run K-means
km <- km.fromscratch(X, 3)
pairs(X, lower.panel = NULL, col = km$cluster)
table(y, km$cluster)

install.packages('mixtools')
library(mixto)

mixmdl = mixtools::normalmixEM(snps$cnlr, k = 2, mu = c(-1,0), ECM = F)


plot_mix_comps = function(x, mu, sigma, lam){
  lam * dnorm(x, mu, sigma)
}
IMPACT.mix.example = data.frame(x = mixmdl$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), 
                 binwidth = 0.05, 
                 colour = "grey55", 
                 fill = "white", 
                 bins = 200,
                 size = 0.1) + 
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                colour = "red", lwd = 1.2) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                colour = "blue", lwd = 1.2) +
  geom_vline(xintercept = mixmdl$mu[which.max(mixmdl$lambda)],
             linetype = 'dashed', 
             size = 0.3, 
             color = 'grey15') +
  scale_y_continuous(expand = c(0.01, 0))
IMPACT.mix.example







##- Nic Socci calls
IM5 = read.csv('06_Nic_Socci_r_001/IM-gbm_IM5_manifest.csv', sep = ',')
IM6 = read.csv('06_Nic_Socci_r_001/IM-gbm_IM6_manifest.csv', sep = ',')
IM7 = read.csv('06_Nic_Socci_r_001/IM-gbm_IM7_manifest.csv', sep = ',')

Nic = rbind(IM5, IM6, IM7)
Nic_tumor = Nic[which(Nic$TYPE == 'Tumor'), ]
write.table(x = Nic_tumor, file = '06_Nic_Socci_r_001/Manifest_Nic.txt', sep = '\t', row.names = F, quote = F)





##-- 
exonic_structure(gene = 'CDKN2A', type = 'Protein coding')
profile()
genecode = rtracklayer::import('~/Documents/GitHub/DryClean_Facets/gencode.v19.chr_patch_hapl_scaff.annotation.gtf.gz')
head(genecode)
selected_gene = genecode %Q% (gene_name == 'CDKN2A')
feature_out = as.data.frame(selected_gene)
feature_out = feature_out[which(feature_out$type == 'protein_coding'), ]
feature_out


ex = feature_out[which(feature_out$type == 'exon'), ]
dim(ex)
ex = ex[!duplicated(ex), ]



EXON1 = chr9 21994138 21994490
EXON2 = chr9 21970901 21971207
EXON3 = chr9 21967751 21968241


countmatrix = readsnpmatrix('01_countmatrices/countsMerged____P-0005701-T01-IM5_P-0005701-N01-IM5.dat.gz/countsMerged____P-0005701-T01-IM5_P-0005701-N01-IM5.dat.gz')
countmatrix = as.data.frame(countmatrix)



out = facetsSuite::run_facets(read_counts = countmatrix,
                              snp_nbhd = 50)

snps = out$snps
snps = snps[which(snps$chrom == 9), ]


a = snps[which(snps$maploc >= 21994138 & snps$maploc <= 21994490), ]
b = snps[which(snps$maploc >= 21970901 & snps$maploc <= 21971207), ]
c =  snps[which(snps$maploc >= 21967751 & snps$maploc <= 21968241), ]

ggplot() +
  geom_raster(data = a, aes(fill = cnlr))





cd = sn[which(sn$maploc >= 21967752 & sn$maploc <= 21995300), ]
cd$s = NA
k = 0
for(i in 1:nrow(cd)){
  cd$s[i] = k
  k = k + 0.1
}

cd$s = cd$s * -1


EXON1 = chr9 21994138 21994490
EXON2 = chr9 21970901 21971207
EXON3 = chr9 21967751 21968241

ggplot(cd) +
  geom_rect(aes(xmin = 21994138, xmax = 21994490, ymin = 0.5, ymax = 1.5)) +
  geom_rect(aes(xmin = 21970901, xmax = 21971207, ymin = 0.5, ymax = 1.5)) +
  geom_rect(aes(xmin = 21967751, xmax = 21968241, ymin = 0.5, ymax = 1.5)) +
  
  geom_point(aes(x = maploc, y = s, color = cnlr), size = 3) +
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





refit = list.files(path = '077_CSF_refit/', recursive = T, pattern = 'gz$', full.names = T)
length(refit)

for(i in 1:length(refit)){
  #print(paste0('mv ', refit[i], ' 08_pileups'))
  system(command = paste0('mv ', refit[i], ' ../08_pileups'))
}


mis = list.files(path = '01_countmatrices/', recursive = T, pattern = 'gz|pileup$', full.names = T)
mis = mis[grep(pattern = 'gz|pileup', mis)]






for file in *.pileup; do 
mv -- "$file" "${file%.pileup}.dat.gz"
done

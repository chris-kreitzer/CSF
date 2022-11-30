##----------------+
## Validation and 
## Justification of the 
## use of Facets and
## chosen parameters
##----------------+

## start: 11/30/2022
## chris-kreitzer

clean()
gc()
.rs.restartR()
setup(working.path = '~/Documents/MSKCC/Subhi/CSF/')
library(data.table)
library(cowplot)

##----------------+
## FACETS vs cBIO (GATK)
## is the usage of FACETS
## justified (overall)?
##----------------+
segmentation = read.csv('Data/Final/71_DMP_cBIO_segmentation.seg', sep = '\t')
sample_pairs = readxl::read_excel('Data/Final/SUBHI SPREADSHEET _USE.xlsx', sheet = 'DMP_CSF_Pairs')
samples_checking = intersect(segmentation$ID, sample_pairs$DMP_CNA_first)
segmentation_cBio = segmentation[which(segmentation$ID %in% samples_checking), ]

##----------------+
## Assess CDKN2A status
## from cBIO calls
##----------------+

#' CDKN2A
CD = data.table(chrom = '9',
                loc.start = 21967751,
                loc.end = 21995300)
key.col = c('chrom', 'loc.start', 'loc.end')

CDKN2A_cbio = data.frame()
for(i in unique(segmentation_cBio$ID)){
  try({
    sample = i
    segs = as.data.table(segmentation_cBio[which(segmentation_cBio$ID == i), ])
    setkey(segs, chrom, loc.start, loc.end)
    overlap = foverlaps(CD, segs, by.x = key.col, by.y = key.col, nomatch = 0)
    patient_seg_mean = overlap$seg.mean
    out = data.frame(sample = sample,
                     chrom = '9',
                     seg_mean = patient_seg_mean)
    CDKN2A_cbio = rbind(CDKN2A_cbio, out)
  })
}


#' EGFR
EG = data.table(chrom = '7',
                loc.start = 55086714,
                loc.end = 55324313)
key.col = c('chrom', 'loc.start', 'loc.end')

EGFR_cbio = data.frame()
for(i in unique(segmentation_cBio$ID)){
  try({
    sample = i
    segs = as.data.table(segmentation_cBio[which(segmentation_cBio$ID == i), ])
    setkey(segs, chrom, loc.start, loc.end)
    overlap = foverlaps(EG, segs, by.x = key.col, by.y = key.col, nomatch = 0)
    patient_seg_mean = overlap$seg.mean
    out = data.frame(sample = sample,
                     chrom = '7',
                     seg_mean = patient_seg_mean)
    EGFR_cbio = rbind(EGFR_cbio, out)
  })
}


#' CDK6
CD6 = data.table(chrom = '7',
                loc.start = 92234235,
                loc.end = 92465908)
key.col = c('chrom', 'loc.start', 'loc.end')

CDK6_cbio = data.frame()
for(i in unique(segmentation_cBio$ID)){
  try({
    sample = i
    segs = as.data.table(segmentation_cBio[which(segmentation_cBio$ID == i), ])
    setkey(segs, chrom, loc.start, loc.end)
    overlap = foverlaps(CD6, segs, by.x = key.col, by.y = key.col, nomatch = 0)
    patient_seg_mean = overlap$seg.mean
    out = data.frame(sample = sample,
                     chrom = '7',
                     seg_mean = patient_seg_mean)
    CDK6_cbio = rbind(CDK6_cbio, out)
  })
}


#' CDK4
CD4 = data.table(chrom = '12',
                 loc.start = 58141510,
                 loc.end = 58149796)
key.col = c('chrom', 'loc.start', 'loc.end')

CDK4_cbio = data.frame()
for(i in unique(segmentation_cBio$ID)){
  try({
    sample = i
    segs = as.data.table(segmentation_cBio[which(segmentation_cBio$ID == i), ])
    setkey(segs, chrom, loc.start, loc.end)
    overlap = foverlaps(CD4, segs, by.x = key.col, by.y = key.col, nomatch = 0)
    patient_seg_mean = overlap$seg.mean
    out = data.frame(sample = sample,
                     chrom = '12',
                     seg_mean = patient_seg_mean)
    CDK4_cbio = rbind(CDK4_cbio, out)
  })
}


#' PTEN
PT = data.table(chrom = '10',
                 loc.start = 89622870,
                 loc.end = 89731687)
key.col = c('chrom', 'loc.start', 'loc.end')

PTEN_cbio = data.frame()
for(i in unique(segmentation_cBio$ID)){
  try({
    sample = i
    segs = as.data.table(segmentation_cBio[which(segmentation_cBio$ID == i), ])
    setkey(segs, chrom, loc.start, loc.end)
    overlap = foverlaps(PT, segs, by.x = key.col, by.y = key.col, nomatch = 0)
    patient_seg_mean = overlap$seg.mean
    out = data.frame(sample = sample,
                     chrom = '10',
                     seg_mean = patient_seg_mean)
    PTEN_cbio = rbind(PTEN_cbio, out)
  })
}



##----------------+
## fetch respective segments 
## from Facets
##----------------+
sample_pairs = readxl::read_excel('Data/Final/SUBHI SPREADSHEET _USE.xlsx', sheet = 'DMP_CSF_Pairs')
sample_pairs = sample_pairs[,c('PatientID', 'DMP_CNA_first', 'CSF_CNA_first')]
sample_pairs = sample_pairs[!is.na(sample_pairs$DMP_CNA_first), ]
GOI = c('CDKN2A', 'EGFR', 'CDK4', 'CDK6', 'PTEN')
folders = list.files(path = '.')
folders = folders[grepl(pattern = 'C-', x = folders)]

DMP_top5 = data.frame()
for(i in 1:nrow(sample_pairs)){
  try({
    DMP = sample_pairs$DMP_CNA_first[i]
    sub_dirs = c(paste0('/Users/chriskreitzer/Documents/MSKCC/Subhi/CSF/', sample_pairs$PatientID[i], '/', DMP, '/'))
    Rfile = list.files(path = sub_dirs, pattern = '.Rdata$', full.names = T)
    
    #' load Facet Fit
    load(file = Rfile)
    
    #' compile the whole FACETS output
    name = DMP
    print(name)
    
    facets_out = list(
      snps = out$jointseg,
      segs = fit$cncf,
      purity = as.numeric(fit$purity),
      ploidy = as.numeric(fit$ploidy),
      dipLogR = out$dipLogR,
      alBalLogR = out$alBalLogR,
      flags = out$flags,
      em_flags = fit$emflags,
      loglik = fit$loglik)
    
    gene_level = facetsSuite::gene_level_changes(facets_output = facets_out, genome = 'hg19')
    
    for(j in GOI){
      print(j)
      if(j %in% unique(gene_level$gene)){
        gene_cnlr = gene_level[which(gene_level$gene == j), 'median_cnlr_seg']
        gene_cnlr = ifelse(is.na(gene_cnlr), NA, gene_cnlr)
        out = data.frame(id = name,
                         gene = j,
                         gene_cnlr = gene_cnlr)
      } else next
      DMP_top5 = rbind(DMP_top5, out)
    }
  })
}

DMP_top5$tag = 'DMP'


##----------------+
## PLOTS
##----------------+
CDKN2A = merge(CDKN2A_cbio, DMP_top5[which(DMP_top5$gene == 'CDKN2A'), ], by.x = 'sample', by.y = 'id', all.x = T)
plotCD2 = ggplot(CDKN2A, aes(x = seg_mean, y = gene_cnlr)) +
  geom_jitter() +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed') +
  scale_y_continuous(limits = c(-5, 1)) +
  scale_x_continuous(limits = c(-5, 1)) +
  theme(aspect.ratio = 1,
        panel.border = element_rect(fill = NA,
                                    linewidth = 1.5, 
                                    color = 'black'),
        axis.text = element_text(size = 10, color = 'black')) +
  labs(x = 'DepthOfCoverage [Cnlr]', y = 'Facets [Cnlr]', title = 'CDKN2A') +
  annotate(geom = 'text', x = -4, y = 0.7, 
           label = paste0('r: ', 
                          round(cor.test(CDKN2A$seg_mean, 
                                         CDKN2A$gene_cnlr, 
                                         method = 'pearson')$estimate[[1]], 3), '\np: ',
                          round(cor.test(CDKN2A$seg_mean, 
                                         CDKN2A$gene_cnlr, 
                                         method = 'pearson')$p.value, 3)))

#' PTEN
PTEN = merge(PTEN_cbio, DMP_top5[which(DMP_top5$gene == 'PTEN'), ], by.x = 'sample', by.y = 'id', all.x = T)
plotPT = ggplot(PTEN, aes(x = seg_mean, y = gene_cnlr)) +
  geom_jitter() +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed') +
  scale_y_continuous(limits = c(-5, 1)) +
  scale_x_continuous(limits = c(-5, 1)) +
  theme(aspect.ratio = 1,
        panel.border = element_rect(fill = NA,
                                    linewidth = 1.5, 
                                    color = 'black'),
        axis.text = element_text(size = 10, color = 'black')) +
  labs(x = 'DepthOfCoverage [Cnlr]', y = 'Facets [Cnlr]', title = 'PTEN') +
  annotate(geom = 'text', x = -4, y = 0.7, 
           label = paste0('r: ', 
                          round(cor.test(PTEN$seg_mean, 
                                         PTEN$gene_cnlr, 
                                         method = 'pearson')$estimate[[1]], 3), '\np: ',
                          round(cor.test(PTEN$seg_mean, 
                                         PTEN$gene_cnlr, 
                                         method = 'pearson')$p.value, 3)))


#' EGFR
EGFR = merge(EGFR_cbio, DMP_top5[which(DMP_top5$gene == 'EGFR'), ], by.x = 'sample', by.y = 'id', all.x = T)
plotEG = ggplot(EGFR, aes(x = seg_mean, y = gene_cnlr)) +
  geom_jitter() +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed') +
  scale_y_continuous(limits = c(-1, 4)) +
  scale_x_continuous(limits = c(-1, 4)) +
  theme(aspect.ratio = 1,
        panel.border = element_rect(fill = NA,
                                    linewidth = 1.5, 
                                    color = 'black'),
        axis.text = element_text(size = 10, color = 'black')) +
  labs(x = 'DepthOfCoverage [Cnlr]', y = 'Facets [Cnlr]', title = 'EGFR') +
  annotate(geom = 'text', x = -0.5, y = 3.7, 
           label = paste0('r: ', 
                          round(cor.test(EGFR$seg_mean, 
                                         EGFR$gene_cnlr, 
                                         method = 'pearson')$estimate[[1]], 3), '\np: ',
                          round(cor.test(EGFR$seg_mean, 
                                         EGFR$gene_cnlr, 
                                         method = 'pearson')$p.value, 3)))

#' CDK4
CDK4 = merge(CDK4_cbio, DMP_top5[which(DMP_top5$gene == 'CDK4'), ], by.x = 'sample', by.y = 'id', all.x = T)
plotCD4 = ggplot(CDK4, aes(x = seg_mean, y = gene_cnlr)) +
  geom_jitter() +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed') +
  scale_y_continuous(limits = c(-1, 4)) +
  scale_x_continuous(limits = c(-1, 4)) +
  theme(aspect.ratio = 1,
        panel.border = element_rect(fill = NA,
                                    linewidth = 1.5, 
                                    color = 'black'),
        axis.text = element_text(size = 10, color = 'black')) +
  labs(x = 'DepthOfCoverage [Cnlr]', y = 'Facets [Cnlr]', title = 'CDK4') +
  annotate(geom = 'text', x = -0.5, y = 3.7, 
           label = paste0('r: ', 
                          round(cor.test(CDK4$seg_mean, 
                                         CDK4$gene_cnlr, 
                                         method = 'pearson')$estimate[[1]], 3), '\np: ',
                          round(cor.test(CDK4$seg_mean, 
                                         CDK4$gene_cnlr, 
                                         method = 'pearson')$p.value, 3)))

#' CDK6
CDK6 = merge(CDK6_cbio, DMP_top5[which(DMP_top5$gene == 'CDK6'), ], by.x = 'sample', by.y = 'id', all.x = T)
plotCD6 = ggplot(CDK4, aes(x = seg_mean, y = gene_cnlr)) +
  geom_jitter() +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed') +
  scale_y_continuous(limits = c(-1, 4)) +
  scale_x_continuous(limits = c(-1, 4)) +
  theme(aspect.ratio = 1,
        panel.border = element_rect(fill = NA,
                                    linewidth = 1.5, 
                                    color = 'black'),
        axis.text = element_text(size = 10, color = 'black')) +
  labs(x = 'DepthOfCoverage [Cnlr]', y = 'Facets [Cnlr]', title = 'CDK6') +
  annotate(geom = 'text', x = -0.5, y = 3.7, 
           label = paste0('r: ', 
                          round(cor.test(CDK6$seg_mean, 
                                         CDK6$gene_cnlr, 
                                         method = 'pearson')$estimate[[1]], 3), '\np: ',
                          round(cor.test(CDK6$seg_mean, 
                                         CDK6$gene_cnlr, 
                                         method = 'pearson')$p.value, 3)))

plot_grid(plotCD2, plotEG, plotCD4, plotCD6, plotPT, ncol = 5, nrow = 1)

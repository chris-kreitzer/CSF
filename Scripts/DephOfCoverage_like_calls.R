##----------------+
## High-level alteration
## Approach alike DepthOfCoverage/IMPACT method
##----------------+

clean()
setwd('~/Documents/MSKCC/11_CSF/')

CSF_samples = read.csv('00_Data/CSF_Samples_to useoncoprint.txt', sep = '\t')
Tumor_samples = read.csv('00_Data/Tumor_Samples_to useoncoprint.txt', sep = '\t')
all_files = list.files('01_countmatrices/', pattern = 'rds$', recursive = T, full.names = T)
GOIs = c('CDKN2A', 'EGFR')

csf_out = data.frame()
for(i in unique(CSF_samples$Sample.ID)){
  print(i)
  file = all_files[grep(pattern = i, x = all_files)]
  file = ifelse(length(file) == 1, file[1], file[grep(pattern = 'adjusted', x = file)])
  datain = readRDS(file = file)
  
  diplogr = datain$dipLogR
  segs = datain$segs
  snps = datain$snps
  
  ##-- Assume that all the CnLR are normally distributed
  Mean = mean(snps$cnlr, na.rm = T)
  Sd = sd(snps$cnlr, na.rm = T)
  
  x = seq(floor(min(snps$cnlr)), ceiling(max(snps$cnlr)), 0.01)
  xp = pnorm(x, mean = Mean, sd = Sd)
  
  ##-- lower threshold 5%
  lower = x[max(which(xp <= 0.1, arr.ind = T))]
  upper = x[min(which(xp >= 0.9, arr.ind = T))]
  
  # plot(x, pnorm(x, mean = Mean, sd = Sd), type = "l",
  #      ylim = c(0, 1), ylab = "P(X < x)", lwd = 2, col = "red")
  # 
  # ##-- lower bound
  # abline(v = x[max(which(xp < 0.05, arr.ind = T))])
  # ##-- upper bound
  # abline(v = x[min(which(xp >= 0.95, arr.ind = T))])
  
  
  ##-- CDKN2A CnLR positions
  cdkn2a = genes_hg19[which(genes_hg19$gene == 'CDKN2A'), ]
  cdkn2a_snps = snps[which(snps$chrom == cdkn2a$chrom &
                             snps$maploc >= cdkn2a$start & 
                             snps$maploc <= cdkn2a$end), ]
  cdkn2a_mean = mean(cdkn2a_snps$cnlr)
  
  ##-- EGFR CnLR positions
  egfr = genes_hg19[which(genes_hg19$gene == 'EGFR'), ]
  egfr_snps = snps[which(snps$chrom == egfr$chrom &
                             snps$maploc >= egfr$start & 
                             snps$maploc <= egfr$end), ]
  egfr_mean = mean(egfr_snps$cnlr)
  
  out = data.frame(id = i,
                   gene = c('CDKN2A', 'EGFR'),
                   mean_value = c(cdkn2a_mean, egfr_mean),
                   n_loci = c(nrow(cdkn2a_snps), nrow(egfr_snps)),
                   lower = lower,
                   upper = upper,
                   call = c(ifelse(cdkn2a_mean <= lower, '-2', 0),
                            ifelse(egfr_mean >= upper, '2', 0)))
  
  csf_out = rbind(csf_out, out)
  rm(segs, snps, datain, file, egfr, cdkn2a, diplogr, lower, upper, x, xp, Mean, Sd)
}
  
 
length(csf_out$id[which(csf_out$gene == 'CDKN2A' & csf_out$call == '-2')]) / length(csf_out$gene[which(csf_out$gene == 'CDKN2A')])
## 24% CDKN2A
length(csf_out$id[which(csf_out$gene == 'EGFR' & csf_out$call == '2')]) / length(csf_out$gene[which(csf_out$gene == 'EGFR')])
## 45% EGFR




##-- TUMOR
tumor_out = data.frame()
for(i in unique(Tumor_samples$Sample.ID)){
  print(i)
  file = all_files[grep(pattern = i, x = all_files)]
  file = ifelse(length(file) == 1, file[1], file[grep(pattern = 'adjusted', x = file)])
  datain = readRDS(file = file)
  
  diplogr = datain$dipLogR
  segs = datain$segs
  
  gene = facetsSuite::gene_level_changes(facets_output = datain, genome = 'hg19', algorithm = 'em')
  gene = gene[which(gene$gene %in% GOIs), c('gene', 'chrom', 'seg', "median_cnlr_seg", 'cn_state')]
  
  segs$adj = ifelse(segs$cnlr.median.clust < 0, segs$cnlr.median.clust - diplogr,
                    segs$cnlr.median.clust - diplogr)
  segs = segs[,c('chrom', 'seg', 'cnlr.median.clust', 'adj')]
  
  gene_new = merge(gene, segs, by = c('chrom', 'seg'), all.x = T)
  
  cnlr_distribution = quantile(segs$adj, probs = c(0.05, 0.25, 0.75, 0.95))
  
  gene_new$call = ifelse(gene_new$adj <= cnlr_distribution[[1]], '-2',
                         ifelse(gene_new$adj >= cnlr_distribution[[4]], '2', 0))
  
  gene_new = gene_new[,c('chrom', 'gene', 'adj', 'call')]
  colnames(gene_new)[3] = 'adjusted_cnlr'
  gene_new$dipLogR = diplogr
  gene_new = gene_new[,c(2,1,3,5,4)]
  gene_new$id = i
  
  tumor_out = rbind(tumor_out, gene_new)
  
  rm(segs, gene, datain, file, gene_new, diplogr, cnlr_distribution)
  
}

length(tumor_out$gene[which(tumor_out$gene == 'CDKN2A' & tumor_out$call == '-2')]) / length(tumor_out$gene[which(tumor_out$gene == 'CDKN2A')])
## 45% CDKN2A deletion
length(tumor_out$gene[which(tumor_out$gene == 'EGFR' & tumor_out$call == '2')]) / length(tumor_out$gene[which(tumor_out$gene == 'EGFR')])
## 49% EGFR amplfications








csf_out = data.frame()
for(i in unique(CSF_samples$Sample.ID)){
  print(i)
  file = all_files[grep(pattern = i, x = all_files)]
  file = ifelse(length(file) == 1, file[1], file[grep(pattern = 'adjusted', x = file)])
  datain = readRDS(file = file)
  
  diplogr = datain$dipLogR
  segs = datain$segs
  snps = datain$snps
  
  gene = facetsSuite::gene_level_changes(facets_output = datain, genome = 'hg19', algorithm = 'em')
  gene = gene[which(gene$gene %in% GOIs), c('gene', 'chrom', 'seg', "median_cnlr_seg", 'cn_state')]
  
  cnlr_distribution = quantile(snps$cnlr, probs = c(0.15, 0.25, 0.75, 0.85))
  segs = segs[,c('chrom', 'seg', 'cnlr.median.clust')]
  
  gene_new = merge(gene, segs, by = c('chrom', 'seg'), all.x = T)
  
  gene_new$call = ifelse(gene_new$median_cnlr_seg <= cnlr_distribution[[1]], '-2',
                         ifelse(gene_new$median_cnlr_seg >= cnlr_distribution[[4]], '2', 0))
  
  gene_new = gene_new[,c('chrom', 'gene', 'median_cnlr_seg', 'call')]
  gene_new$lower_cutoff = cnlr_distribution[[1]]
  gene_new$upper_cutoff = cnlr_distribution[[4]]
  gene_new$id = i
  
  csf_out = rbind(csf_out, gene_new)
  
  rm(segs, gene, datain, file, gene_new, diplogr, cnlr_distribution)
  
}

length(csf_out$gene[which(csf_out$gene == 'CDKN2A' & csf_out$call == '-2')]) / length(csf_out$gene[which(csf_out$gene == 'CDKN2A')])
## 15% CDKN2A
length(csf_out$gene[which(csf_out$gene == 'EGFR' & csf_out$call == '2')]) / length(csf_out$gene[which(csf_out$gene == 'EGFR')])
## 21% EGFR

tumor_out = data.frame()
for(i in unique(Tumor_samples$Sample.ID)){
  print(i)
  file = all_files[grep(pattern = i, x = all_files)]
  file = ifelse(length(file) == 1, file[1], file[grep(pattern = 'adjusted', x = file)])
  datain = readRDS(file = file)
  
  diplogr = datain$dipLogR
  segs = datain$segs
  snps = datain$snps
  
  gene = facetsSuite::gene_level_changes(facets_output = datain, genome = 'hg19', algorithm = 'em')
  gene = gene[which(gene$gene %in% GOIs), c('gene', 'chrom', 'seg', "median_cnlr_seg", 'cn_state')]
  
  cnlr_distribution = quantile(snps$cnlr, probs = c(0.05, 0.25, 0.75, 0.95))
  segs = segs[,c('chrom', 'seg', 'cnlr.median.clust')]
  
  gene_new = merge(gene, segs, by = c('chrom', 'seg'), all.x = T)
  
  gene_new$call = ifelse(gene_new$median_cnlr_seg <= cnlr_distribution[[1]], '-2',
                         ifelse(gene_new$median_cnlr_seg >= cnlr_distribution[[4]], '2', 0))
  
  gene_new = gene_new[,c('chrom', 'gene', 'median_cnlr_seg', 'call')]
  gene_new$lower_cutoff = cnlr_distribution[[1]]
  gene_new$upper_cutoff = cnlr_distribution[[4]]
  gene_new$id = i
  
  tumor_out = rbind(tumor_out, gene_new)
  
  rm(segs, gene, datain, file, gene_new, diplogr, cnlr_distribution)
  
}

length(tumor_out$gene[which(tumor_out$gene == 'CDKN2A' & tumor_out$call == '-2')]) / length(tumor_out$gene[which(tumor_out$gene == 'CDKN2A')])
## 34% CDKN2A deletion
length(tumor_out$gene[which(tumor_out$gene == 'EGFR' & tumor_out$call == '2')]) / length(tumor_out$gene[which(tumor_out$gene == 'EGFR')])
## 30% EGFR amplfications









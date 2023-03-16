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
  
  gene = facetsSuite::gene_level_changes(facets_output = datain, genome = 'hg19', algorithm = 'em')
  gene = gene[which(gene$gene %in% GOIs), c('gene', 'chrom', 'seg', "median_cnlr_seg", 'cn_state')]
  
  segs$adj = ifelse(segs$cnlr.median.clust < 0, segs$cnlr.median.clust - diplogr,
                    segs$cnlr.median.clust - diplogr)
  segs = segs[,c('chrom', 'seg', 'cnlr.median.clust', 'adj')]
  
  gene_new = merge(gene, segs, by = c('chrom', 'seg'), all.x = T)
  
  cnlr_distribution = quantile(segs$adj, probs = c(0.15, 0.25, 0.75, 0.85))
  
  gene_new$call = ifelse(gene_new$adj <= cnlr_distribution[[1]], '-2',
                         ifelse(gene_new$adj >= cnlr_distribution[[4]], '2', 0))
  
  gene_new = gene_new[,c('chrom', 'gene', 'adj', 'call')]
  colnames(gene_new)[3] = 'adjusted_cnlr'
  gene_new$dipLogR = diplogr
  gene_new = gene_new[,c(2,1,3,5,4)]
  gene_new$id = i
  
  csf_out = rbind(csf_out, gene_new)
  
  rm(segs, gene, datain, file, gene_new, diplogr, cnlr_distribution)
  
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

















GBM/07289_O/1.1.2/20220329_19_53_495980/bam/s_C_RT2H2R_N901_dZ_IM6.rg.md.abra.printreads.bam

























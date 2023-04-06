##----------------+
## High-level alteration
## Approach like DepthOfCoverage/IMPACT method
##----------------+


## start: 03/28/2023
## revision: 04/05/2023
## revision: 04/06/2023
## chris-kreitzer


clean()
setwd('~/Documents/MSKCC/11_CSF/')

database = readxl::read_excel('00_Data/Chris_Subhi_Data_April_5_2023.xlsx')
files = list.files('01_countmatrices/', pattern = 'gz|pileup', recursive = T, full.names = T)
GOIs = c('CDKN2A','CDK4','CDK6','PTEN','EGFR','PDGFRA','KIT','KDR','MET','MDM2','MDM4','RB1','NF1','TP53','FGF4', 'FGF19')

paths = data.frame()
for(i in 1:nrow(database)){
  if(any(grepl(pattern = database$Sample.ID[i], x = files))){
    path = files[grep(pattern = database$Sample.ID[i], x = files)]
    path = ifelse(length(path) != 1, path[grep(pattern = 'gz|pileup', x = path)], path[1])
    out = data.frame(id = database$Sample.ID[i],
                     countmatrix = path)
  } else {
    print(database$Sample.ID[i])
  }
  paths = rbind(paths, out)
}


##-- determine the accessibility of count-matrices
err.thresh = 10
del.thresh = 10

all_out = data.frame()
for(i in 1:nrow(paths)){
  print(i)
  try({
    if(any(grepl(pattern = 'gz', x = paths$countmatrix[i]))){
      pileup = facetsSuite::read_snp_matrix(input_file = paths$countmatrix[i])
      countmatrix = as.data.frame(pileup[,c(1,2,3,5,4,6)])
    } else {
      pileup = read.csv(file = paths$countmatrix[i], sep = ',')
      ii = which(pileup$File1E <= err.thresh & pileup$File1D <= del.thresh & pileup$File2E <= err.thresh & pileup$File2D <= del.thresh)
      rcmat = pileup[ii, 1:2]
      rcmat$NOR.DP = pileup$File1R[ii] + pileup$File1A[ii]
      rcmat$NOR.RD = pileup$File1R[ii]
      rcmat$TUM.DP = pileup$File2R[ii] + pileup$File2A[ii]
      rcmat$TUM.RD = pileup$File2R[ii]
      countmatrix = as.data.frame(rcmat)
    }
    
    ##-- run Facets and fetch SNP database
    out = facetsSuite::run_facets(read_counts = countmatrix, 
                                  cval = 100,
                                  dipLogR = NULL,
                                  snp_nbhd = 250,
                                  seed = 100, 
                                  genome = 'hg19', 
                                  ndepth = 20)
    ##-- investigate pileup files
    snps = out$snps
    pdf(file = paste0(paths$countmatrix[i], '_Allelic_Distribution.pdf'), width = 4.5, height = 4.5)
    plot(density(snps$vafT[which(snps$het == 1)]), main = 'allele freq. distribution het. SNPs')
    dev.off()
    
    if(length(snps$maploc[which(snps$het == 1)]) <= 1000 | 
       dim(snps)[1] < 5000){
      flag_pileup = T
    } else {
      flag_pileup = F
    }
    
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
    for(gene in unique(GOIs)){
      position = snps[which(snps$chrom == genes_hg19$chrom[which(genes_hg19$gene == gene)] &
                              snps$maploc >= genes_hg19$start[which(genes_hg19$gene == gene)] &
                              snps$maploc <= genes_hg19$end[which(genes_hg19$gene == gene)]), ]
      gene_cnlr_mean = mean(position$cnlr)
      gene_call = ifelse(gene_cnlr_mean >= upper, 2, 
                         ifelse(gene_cnlr_mean <= lower, -2, 0))
      hugo = gene
      
      gene_out = data.frame(id = paths$id[i],
                            gene = hugo,
                            call = gene_call,
                            upper = upper,
                            lower = lower,
                            n_hets = length(snps$maploc[which(snps$het == 1)]),
                            flag_pileup = flag_pileup)
      
      all_out = rbind(all_out, gene_out)
    }
    
    rm(out, countmatrix, lower, upper, position, gene, gene_cnlr_mean, gene_out, Mean, Sd, x, xp)
  })
}

write.table(x = all_out, file = '00_Data/MSK_like_CNA_calls.txt', sep = '\t', row.names = F, quote = F)



paths$id[81]



out = facetsSuite::read_snp_matrix('01_countmatrices/countsMerged____P-0000651-T01-IM3_P-0000651-N01-IM3.dat.gz/countsMerged____P-0000651-T01-IM3_P-0000651-N01-IM3.dat.gz')
out = as.data.frame(out[,c(1,2,3,5,4,6)])
ss = facetsSuite::run_facets(read_counts = out)
snps = ss$snps
pdf(file = paste0(paths$countmatrix[1], '_allelic_Distribution.pdf'), width = 4, height = 4)
plot(density(snps$vafT[which(snps$het == 1)]))
dev.off()










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
  
 
length(unique(csf_out$id[which(csf_out$gene == 'CDKN2A' & csf_out$call == '-2')])) / length(csf_out$gene[which(csf_out$gene == 'CDKN2A')])
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









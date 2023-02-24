##----------------+
## Reduced fits
## - check the samples that failed the full
## Facets approach due to purity, and/or icn filter
## - rework them, and get absolute numbers
## based on cnlr_median_clust
##----------------+


clean()
gc()
setwd('~/Documents/MSKCC/11_CSF/01_countmatrices/')
load('../test_facets_output.rda')
source('~/Documents/MSKCC/11_CSF/02_Scripts/cf_Plot.R')
library(facetsSuite)
library(facets)
library(patchwork)


GOIs = c('CDKN2A','CDK4','CDK6','PTEN','EGFR','RET',
         'GATA3','PDGFRA','KIT','KDR','MET','MDM2',
         'MDM4','RB1','TP53')
copynumberstates = facetsSuite:::copy_number_states

database = readxl::read_excel('../00_Data/Database_Final+Feb23_ck_2023.xlsx', sheet = 1)


##----------------+
## Work on FULL FACETS samples
## - purity, ploidy, 
## - FGA, oncogenic genes
## - and 10p and 7q
## - verify with cBIO Portal for DMP samples
##----------------+
full_facets = database[which(database$Fit == 'full'), ]
full_facets$id = basename(full_facets$FacetsCountfile)
full_facets = as.data.frame(full_facets)

alteration_summary = data.frame()
gene_calls_out = data.frame()

for(i in 1:nrow(full_facets)){
  print(i)
  files = list.files(path = full_facets$id[i], full.names = T, recursive = T)
  fit = readRDS(file = files[grep(pattern = '.rds$', x = files)])
  gene = read.csv(file = files[grep(pattern = '_gene_level.txt$', x = files)], sep = '\t')
  qc = facets_fit_qc(facets_output = fit)
  arm_change = facetsSuite::arm_level_changes(segs = fit$segs,
                                              ploidy = fit$ploidy,
                                              genome = 'hg19',
                                              algorithm = 'em')
  
  ##-- General measures
  purity = qc$purity
  ploidy = qc$ploidy
  wgd = qc$wgd
  fga = qc$fga
  
  
  ##-- assign numeric call for arms
  arm_10p = arm_change$full_output[which(arm_change$full_output$arm == '10p'), c('arm', 'cn_state')]
  if(nrow(arm_10p) == 0){
    arm_10p = NA
  } else {
    arm_10p = unique(copynumberstates$numeric_call[which(copynumberstates$call == arm_10p$cn_state)])
    arm_10p = ifelse(length(arm_10p) == 2 & wgd, min(arm_10p), max(arm_10p))
  }
  
  arm_10q = arm_change$full_output[which(arm_change$full_output$arm == '10q'), c('arm', 'cn_state')]
  if(nrow(arm_10q) == 0){
    arm_10q = NA
  } else {
    arm_10q = unique(copynumberstates$numeric_call[which(copynumberstates$call == arm_10q$cn_state)])
    arm_10q = ifelse(length(arm_10q) == 2 & wgd, min(arm_10q), max(arm_10q))
  }
  
  arm_7p = arm_change$full_output[which(arm_change$full_output$arm == '7p'), c('arm', 'cn_state')]
  if(nrow(arm_7p) == 0){
    arm_7p = NA
  } else {
    arm_7p = unique(copynumberstates$numeric_call[which(copynumberstates$call == arm_7p$cn_state)])
    arm_7p = ifelse(length(arm_7p) == 2 & wgd, min(arm_7p), max(arm_7p))
  }
  
  arm_7q = arm_change$full_output[which(arm_change$full_output$arm == '7q'), c('arm', 'cn_state')]
  if(nrow(arm_7q) == 0){
    arm_7q = NA
  } else {
    arm_7q = unique(copynumberstates$numeric_call[which(copynumberstates$call == arm_7q$cn_state)])
    arm_7q = ifelse(length(arm_7q) == 2 & wgd, min(arm_7q), max(arm_7q))
  }
  
  
  ##-- Gene Level alterations
  
  onco_cna = gene[which(gene$gene %in% GOIs), c('gene', 'tcn.em', 'mcn', 'lcn.em', 'cn_state')]
  onco_cna$string = paste(wgd, onco_cna$tcn.em, onco_cna$mcn, onco_cna$lcn.em, sep = ':')
  onco_cna = merge(onco_cna, copynumberstates[,c('numeric_call', 'map_string')],
                   by.x = 'string', by.y = 'map_string', all.x = T)
  onco_cna$numeric_call = ifelse(is.na(onco_cna$numeric_call) & onco_cna$tcn.em > 5, 2, onco_cna$numeric_call)
  onco_cna$id = full_facets$id[i]
  
  out = data.frame(id = full_facets$id[i],
                   arm_7p = arm_7p,
                   arm_7q = arm_7q,
                   arm_10p = arm_10p,
                   arm_10q = arm_10q,
                   purity = purity,
                   ploidy = ploidy,
                   wgd = wgd,
                   fga = fga,
                   oncogenic_CNA = paste(onco_cna$gene[which(onco_cna$numeric_call %in% c(-2, -1, 1, 2))], collapse = ','))
  
  alteration_summary = rbind(alteration_summary, out)
  gene_calls_out = rbind(gene_calls_out, onco_cna)
  rm(gene,qc, arm_10p, arm_10q, arm_7p, arm_7q, arm_change, fit, purity, ploidy, fga)

}

write.table(x = alteration_summary, file = '../00_Data/FULL_FACETS_sampleSummary.txt', sep = '\t', row.names = F, quote = F)
write.table(x = gene_calls_out, file = '../00_Data/FULL_FACETS_geneLevelSummary.txt', sep = '\t', row.names = F, quote = F)




##----------------+
## REDUCED FACETS
## - fetch GOI alteration signals
## - no purity/ploidy/fga estimates
##----------------+
reduced_facets = database[which(database$Fit == 'reduced'), ]
reduced_facets$id = basename(reduced_facets$FacetsCountfile)
reduced_facets = as.data.frame(reduced_facets)


##-- New fits and CnLR of selected genes
err.thresh = 10
del.thresh = 10

number = 79
basename(reduced_facets$id[number])
samplePath = paste0(reduced_facets$id[number], '/', 
                    list.files(path = reduced_facets$id[number], pattern = 'dat.gz$', full.names = F)) 
sampleid = 's_PN_FROZEN_07289_H_cas.rg.md.abra.printreads__s_C_006878_S001_d01'


countmatrix = facetsSuite::read_snp_matrix(input_file = samplePath)
pileup = read.csv(file = samplePath, sep = ',')
ii = which(pileup$File1E <= err.thresh & pileup$File1D <= del.thresh & pileup$File2E <= err.thresh & pileup$File2D <= del.thresh)
rcmat = pileup[ii, 1:2]
rcmat$NOR.DP = pileup$File1R[ii] + pileup$File1A[ii]
rcmat$NOR.RD = pileup$File1R[ii]
rcmat$TUM.DP = pileup$File2R[ii] + pileup$File2A[ii]
rcmat$TUM.RD = pileup$File2R[ii]
countmatrix = rcmat

countmatrix = countmatrix[,c(1,2,3,5,4,6)]
fit = facetsSuite::run_facets(read_counts = countmatrix, 
                              cval = 40,
                              dipLogR = 0.39,
                              snp_nbhd = 250,
                              seed = 100)
i = facetsSuite::cnlr_plot(fit, return_object = T)
ii = facetsSuite::valor_plot(fit, return_object = T)
iii = facetsSuite::icn_plot(fit, return_object = T)
iv = cf_plot(facets_data = fit, method = 'em', genome = 'hg19')
i/ii/iii/iv +
  plot_layout(heights = c(1,1,1,0.25))


diplogr = fit$dipLogR
segs = fit$segs


gene = facetsSuite::gene_level_changes(facets_output = fit, genome = 'hg19', algorithm = 'em')
gene = gene[which(gene$gene %in% GOIs), c('gene', 'chrom', 'seg', "median_cnlr_seg", 'cn_state')]
segs$adj = ifelse(segs$cnlr.median.clust < 0, segs$cnlr.median.clust - diplogr,
                  segs$cnlr.median.clust - diplogr)
segs = segs[,c('chrom', 'seg', 'cnlr.median.clust', 'adj')]

gene_new = merge(gene, segs, by = c('chrom', 'seg'), all.x = T)

cnlr_distribution = quantile(segs$adj, probs = c(0.1, 0.25, 0.75, 0.9))

gene_new$call = ifelse(gene_new$adj <= cnlr_distribution[[1]], '-2',
                   ifelse(gene_new$adj >= cnlr_distribution[[4]], '2', 0))

gene_new = gene_new[,c('chrom', 'gene', 'adj', 'call')]
colnames(gene_new)[3] = 'adjusted_cnlr'
gene_new$dipLogR = diplogr
gene_new = gene_new[,c(2,1,3,5,4)]
gene_new


dir.create(path = paste0(basename(samplePath), '/', 'adjusted/'))
write.table(x = gene_new, file = paste0(basename(samplePath), '/', 'adjusted/', sampleid, '_gene_level.txt'), 
            sep = '\t', quote = F, row.names = F)
saveRDS(object = fit, file = paste0(basename(samplePath), '/', 'adjusted/', sampleid, '.rds'))

rm(fit, countmatrix, i, ii, iii, iv, diplogr, segs, gene, gene_new, cnlr_distribution, sampleid, samplePath)




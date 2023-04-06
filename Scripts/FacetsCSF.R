##----------------+
## FACETS fits for Tumor/CSF samples
## Check all_files - location + saving
##----------------+
##
## start: 02/11/2023
## revision: 02/22/2023
## revision: 02/24/2023
## chris-kreitzer



setwd('~/Documents/MSKCC/11_CSF/01_countmatrices/')
source('~/Documents/MSKCC/11_CSF/02_Scripts/cf_Plot.R')
library(patchwork)
library(facets)
library(pctGCdata)
library(facetsSuite)

GOIs = c('CDKN2A','CDK4','CDK6','PTEN','EGFR','PDGFRA','KIT','KDR','MET','MDM2','MDM4','RB1','NF1','TP53','FGF4', 'FGF19')
all_files = list.files('../03_missing/', full.names = T, recursive = T)
err.thresh = 10
del.thresh = 10


##-- START
basename(all_files[1])
samplePath = all_files[1]
sampleid = 'countsMerged____P-0063687-T01-IM7_P-0063687-N01-IM7'


countmatrix = facetsSuite::read_snp_matrix(input_file = samplePath)
# pileup = read.csv(file = samplePath, sep = ',')
# ii = which(pileup$File1E <= err.thresh & pileup$File1D <= del.thresh & pileup$File2E <= err.thresh & pileup$File2D <= del.thresh)
# rcmat = pileup[ii, 1:2]
# rcmat$NOR.DP = pileup$File1R[ii] + pileup$File1A[ii]
# rcmat$NOR.RD = pileup$File1R[ii]
# rcmat$TUM.DP = pileup$File2R[ii] + pileup$File2A[ii]
# rcmat$TUM.RD = pileup$File2R[ii]
# countmatrix = rcmat

countmatrix = countmatrix[,c(1,2,3,5,4,6)]
out = facetsSuite::run_facets(read_counts = countmatrix, 
                              cval = 100,
                              dipLogR = NULL,
                              snp_nbhd = 250,
                              seed = 100)
i = facetsSuite::cnlr_plot(out, return_object = T)
ii = facetsSuite::valor_plot(out, return_object = T)
iii = facetsSuite::icn_plot(out, return_object = T)
iv = cf_plot(facets_data = out, method = 'em', genome = 'hg19')
qc = facets_fit_qc(out)
qc$facets_qc
i/ii/iii/iv +
  plot_layout(heights = c(1,1,1,0.25))
qc$ploidy
qc$purity
qc$wgd
out$dipLogR


##-- FINAL
fit = facetsSuite::run_facets(read_counts = countmatrix, 
                              cval = 50,
                              dipLogR = -0.058,
                              snp_nbhd = 250,
                              seed = 100)
i = facetsSuite::cnlr_plot(fit, return_object = T)
ii = facetsSuite::valor_plot(fit, return_object = T)
iii = facetsSuite::icn_plot(fit, return_object = T)
iv = cf_plot(facets_data = fit, method = 'em', genome = 'hg19')
qc = facets_fit_qc(fit)
qc$facets_qc
i/ii/iii/iv +
  plot_layout(heights = c(1,1,1,0.25))
qc$ploidy
qc$purity
qc$wgd
fit$dipLogR


##--
gene = facetsSuite::gene_level_changes(facets_output = fit)
gene[which(gene$gene %in% GOIs), c('gene', 'tcn.em', 'lcn.em', 'cn_state', 'filter')]
qc = as.data.frame(qc)
cncf = fit$segs
png_a = i/ii/iii/iv +
  plot_layout(heights = c(1,1,1,0.25))
IGV = facetsSuite::format_igv_seg(facets_output = fit, sample_id = sampleid, normalize = T)


# dir.create(path = paste0('../01_countmatrices/', sub('^(.*[\\/])',"", samplePath)))
# dir_new = paste0('../01_countmatrices/', sub('^(.*[\\/])',"", samplePath))
# system(command = paste0('mv ', samplePath, ' ', dir_new))


write.table(x = gene, file = paste0(sub('^(.*[\\/])',"", samplePath), '/', sampleid, '_gene_level.txt'), 
            sep = '\t', quote = F, row.names = F)
write.table(x = cncf, file = paste0(sub('^(.*[\\/])',"", samplePath), '/', sampleid, '_cncf.txt'), 
            sep = '\t', quote = F, row.names = F)
write.table(x = qc, file = paste0(sub('^(.*[\\/])',"", samplePath), '/', sampleid, '_qc.txt'), 
            sep = '\t', quote = F, row.names = F)
write.table(x = IGV, file = paste0(sub('^(.*[\\/])',"", samplePath), '/', sampleid, '_adjusted.seg'), 
            sep = '\t', quote = F, row.names = F)

saveRDS(object = fit, file = paste0(sub('^(.*[\\/])',"", samplePath), '/', sampleid, '.rds'))
ggsave(filename = paste0(sub('^(.*[\\/])',"", samplePath), '/', sampleid, '.png'), 
       plot = png_a, device = 'png')

rm(out, fit, qc, cncf, gene, IGV, png_a, countmatrix, pileup, rcmat)


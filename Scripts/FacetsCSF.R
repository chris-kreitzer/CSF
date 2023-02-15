setwd('~/Documents/MSKCC/11_CSF/01_countmatrices/')
library(patchwork)
library(facets)
library(pctGCdata)
library(facetsSuite)

facets = read.csv('~/Documents/MSKCC/05_IMPACT40K/Data/Signed_out/Facets_annotated.cohort.txt', sep = '\t')
GOIs = c('CDKN2A','CDK4','CDK6','PTEN','EGFR','PDGFRA','KIT','KDR','MET','MDM2','MDM4','RB1','NF1','TP53','FGF4', 'FGF19')
all_files = list.files('.', full.names = T, recursive = T)


##-- START
basename(all_files[27])
samplePath = all_files[27]
sampleid = 'P-0005361-T01-IM5_P-0005361-N01-IM5'

countmatrix = facetsSuite::read_snp_matrix(input_file = samplePath)
countmatrix = countmatrix[,c(1,2,3,5,4,6)]
out = facetsSuite::run_facets(read_counts = countmatrix, 
                              cval = 100,
                              dipLogR = -0.03,
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
                              cval = 90,
                              dipLogR = -0.03,
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

rm(out, fit, qc, cncf, gene, IGV, png_a)

## Automated pipeline for CNA detection
## by using Facets
## 
## start: 09/11/2022
## revision: 09/28/2022
## revision: 11/25/2022
## revision: 04/15/2023
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
source('~/Documents/GitHub/CSF/Scripts/UtilityFunctions.R')
source('~/Documents/GitHub/CSF/Scripts/FacetsPlot.R')
source('~/Documents/GitHub/CSF/Scripts/cf_Plot.R')
source('~/Documents/GitHub/CSF/Scripts/ReadSnpMatrix.R')


database = readxl::read_excel('00_Data/Database_Chris_Sub_April12.xlsx')
csf = database[which(database$TYPE == 'CSF'), ]
csf = csf[which(csf$CSF_STATUS != 'Test Failure'), ]
files = list.files(path = '08_pileups/', full.names = T)


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
##----------------+
number = 155
sample = csf$Sample.ID[number]

countmatrix = readsnpmatrix(path = files[grep(pattern = sample, x = files)])
countmatrix = as.data.frame(countmatrix)

snps = facetsSuite::run_facets(read_counts = countmatrix,
                               cval = 100,
                               dipLogR = NULL,
                               snp_nbhd = 250,
                               seed = 100, 
                               genome = 'hg19', 
                               ndepth = 20)
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



##----------------+
## Second run:
## Ploidy, Purity, dipLogR determination
## (exclusively purity runs);
##----------------+
cval = 100
seed = 100
min_het = 15
genome = 'hg19'
snp_nbhd = 250

out = facetsSuite::run_facets(read_counts = countmatrix,
                               cval = cval,
                               dipLogR = NULL,
                               snp_nbhd = snp_nbhd,
                               seed = seed, 
                               genome = 'hg19', 
                               ndepth = 30)
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

sample_summary = data.frame(id = sample,
                            CNA_fit = 'pass',
                            Flag = '')

write.table(x = sample_summary, file = paste0('07_CSF_refit/', sample, '/', sample, '_summary.txt'), sep = '\t', row.names = F, quote = F)
saveRDS(object = out, file = paste0('07_CSF_refit/', sample, '/', sample, '_second_pass.rds'))
out$dipLogR
rm(i, ii, iii, iv, pass2, sample_summary, out, cval, min_het, het_snps, snps)




##----------------+
## Third run:
## Gene-Level alterations if all the steps before
## are satisfied
## in-depth analyzes of CKDN2A, CDK6, and EGFR
##----------------+
cval = 50
seed = 100
min_het = 15
genome = 'hg19'
snp_nbhd = 100
diplogr = 0.1281249


fit = facetsSuite::run_facets(read_counts = countmatrix,
                              cval = 40,
                              dipLogR = diplogr,
                              snp_nbhd = snp_nbhd,
                              seed = seed, 
                              genome = 'hg19', 
                              ndepth = 30)

i = cnlr_plot(facets_data = fit, genome = 'hg19')
ii = valor_plot(facets_data = fit, genome = 'hg19')
iii = icn_plot(facets_data = fit, genome = 'hg19')
iv = cf_plot(facets_data = fit, genome = 'hg19')

i / ii / iii / iv + plot_layout(heights = c(1,1,0.5,0.25))

qc = facets_fit_qc(facets_output = fit)
qc


## Genes:
cdkn2a_first = fit$snps[which(fit$snps$chrom == 9 &
                          fit$snps$maploc >= 21967751 &
                          fit$snps$maploc <= 21995300), ]
cdkn2a_first
cdkn2a_second = fit$segs[which(fit$segs$seg == unique(cdkn2a_first$seg)), ]
cdkn2a_second


egfr_first = fit$snps[which(fit$snps$chrom == 7 &
                                fit$snps$maploc >= 55086714 &
                                fit$snps$maploc <= 55324313), ]
egfr_second = fit$segs[which(fit$segs$seg == unique(egfr_first$seg)), ]
egfr_second


genes_all = facetsSuite::gene_level_changes(facets_output = fit, genome = 'hg19')
genes_all[which(genes_all$gene == 'EGFR'), c('chrom', 'median_cnlr_seg', 'tcn.em', 'lcn.em', 'cn_state', 'filter')]
genes_all[which(genes_all$gene == 'CDKN2A'), c('chrom', 'median_cnlr_seg', 'tcn.em', 'lcn.em', 'cn_state', 'filter')]
genes_all[which(genes_all$gene == 'CDK6'), c('chrom', 'median_cnlr_seg', 'tcn.em', 'lcn.em', 'cn_state', 'filter')]

write.table(x = genes_all, file = paste0('07_CSF_refit/', sample, '/', sample, '_gene_level_alteration.txt'), sep = '\t', row.names = F)
saveRDS(object = fit, file = paste0('07_CSF_refit/', sample, '/', sample, '_third_pass.rds'))

rm(cdkn2a_first, cdkn2a_second, egfr_first, egfr_second, fit, qc, norm_density, countmatrix, i, ii, iii, iv)


## out























#' ## Parameters: (exclusively purity runs); not interested in gene_level alterations
#' cval = 100
#' seed = 100
#' min_het = 15
#' genome = 'hg19'
#' 
#' 
#' ##-----------------
#' ## First RUN
#' ##-----------------
#' parameter_table = data.frame(tumor_sample = snp_pileup$pileup_file[which(snp_pileup$Patient_ID == ID)],
#'                              name = basename(snp_pileup$sample[which(snp_pileup$Patient_ID == ID)]),
#'                              dipLogR = NA)
#' 
#' if(any(grepl(pattern = '_N90', parameter_table$name) | grepl(pattern = 'POOLED', parameter_table$name) | 
#'        grepl(pattern = '_FROZEN', parameter_table$name) | grepl(pattern = '_Frozen', parameter_table$name) | 
#'        grepl(pattern = '_NCAS', parameter_table$name) |
#'        grepl(pattern = '_FFPE', parameter_table$name))){
#'   parameter_table = parameter_table[!grepl(pattern = '_N90', parameter_table$name), ]
#'   parameter_table = parameter_table[!grepl(pattern = 'POOLED', parameter_table$name), ]
#'   parameter_table = parameter_table[!grepl(pattern = '_FROZEN', parameter_table$name), ]
#'   parameter_table = parameter_table[!grepl(pattern = '_Frozen', parameter_table$name), ]
#'   parameter_table = parameter_table[!grepl(pattern = '_NCAS', parameter_table$name), ]
#'   parameter_table = parameter_table[!grepl(pattern = '_FFPE', parameter_table$name), ]
#'   parameter_table$tumor_sample = seq(2, nrow(parameter_table) + 1, 1)
#' }
#' 
#' 
#' dev.off()
#' gene_level_out = data.frame()
#' facets_plots = list()
#' for(tumor_sample in 1:nrow(parameter_table)){
#'   countMatrix = multi_readSnpMatrix(filename = countMatrix_path, 
#'                                     tumor_sample = parameter_table$tumor_sample[tumor_sample])
#'   fit = facetsSuite::run_facets(read_counts = countMatrix,
#'                                   cval = cval,
#'                                   min_nhet = min_het,
#'                                   genome = genome,
#'                                   seed = seed)
#'   
#'   print(paste0(parameter_table$name[tumor_sample], ': facets_qc ', facets_fit_qc(fit)$facets_qc))
#'   filters = c('waterfall_filter_pass', 'hyper_seg_filter_pass',
#'               'em_cncf_icn_discord_filter_pass', 'contamination_filter_pass')
#'   facets_qc = facets_fit_qc(fit)
#'   
#'   
#'   if(all(unlist(facets_qc[filters]))){
#'     qc_5parameters = TRUE
#'     print(paste0('min qc :', qc_5parameters))
#'   } else {
#'     qc_5parameters = FALSE
#'     print(paste0('min qc :', qc_5parameters))
#'   }
#'   
#'   i = facetsSuite::cnlr_plot(fit, return_object = T)
#'   ii = facetsSuite::valor_plot(fit, return_object = T)
#'   iii = facetsSuite::icn_plot(fit, return_object = T)
#'   iv = facetsSuite::cf_plot(fit, return_object = T)
#'   print(i / ii/ iii/ iv)
#'   all_plots = i / ii / iii / iv
#'   
#'   #' gene level changes
#'   gene_level = facetsSuite::gene_level_changes(fit,
#'                                                genome = 'hg19')
#'   dipLogR_original = fit$dipLogR
#'   print(paste0('Original dipLogR: ', dipLogR_original))
#'   purity = fit$purity
#'   
#'   #' run on clinically significant genes
#'   for(gene in c('CDKN2A',
#'                 'CDK4',
#'                 'CDK6',
#'                 'PTEN',
#'                 'EGFR',
#'                 'PDGFRA',
#'                 'KIT',
#'                 'KDR',
#'                 'MET',
#'                 'MDM2',
#'                 'MDM4',
#'                 'RB1',
#'                 'NF1',
#'                 'TP53')){
#'       try({
#'         #' Facets CnLR estimates
#'         if(facets_qc$facets_qc & qc_5parameters){
#'           if(gene %in% gene_level$gene){
#'             CnLR = gene_level[which(gene_level$gene == gene), 'median_cnlr_seg']
#'             CnLR = ifelse(is.nan(CnLR) | is.na(CnLR) | length(CnLR) == 0, NA, CnLR - dipLogR_original)
#'             cf.em = gene_level[which(gene_level$gene == gene), 'cf.em']
#'             cf.em = ifelse(is.nan(cf.em) | is.na(cf.em) | length(cf.em) == 0, NA, cf.em)
#'             gene_cn = gene_level[which(gene_level$gene == gene), 'cn_state']
#'             clonality = ifelse(cf.em >= (purity * 0.8), 'clonal', 'subclonal')
#'             pass = gene_level[which(gene_level$gene == gene), 'filter']
#'             
#'           } else {
#'             CnLR = NA
#'             cf.em = NA
#'             gene_cn = NA
#'             clonality = NA
#'             pass = NA
#'           }
#'           
#'         } else if (facets_qc$facets_qc | qc_5parameters) {
#'           if(gene %in% gene_level$gene){
#'             gene_chrom = as.integer(genes_hg19[which(genes_hg19$gene == gene), 'chrom'])
#'             gene_start = as.integer(genes_hg19[which(genes_hg19$gene == gene), 'start'])
#'             gene_end = as.integer(genes_hg19[which(genes_hg19$gene == gene), 'end'])
#'             gene_snps = fit$snps[which(fit$snps$chrom == gene_chrom & 
#'                                          fit$snps$maploc >= gene_start & 
#'                                          fit$snps$maploc <= gene_end), 'cnlr']
#'             gene_snps = as.numeric(gene_snps)
#'             CnLR = ifelse(is.nan(gene_snps) | is.na(gene_snps) | length(gene_snps) == 0, NA, mean(gene_snps, na.rm = T))
#'             cf.em = gene_level[which(gene_level$gene == gene), 'cf.em']
#'             cf.em = ifelse(is.nan(cf.em) | is.na(cf.em) | length(cf.em) == 0, NA, cf.em)
#'             gene_cn = gene_level[which(gene_level$gene == gene), 'cn_state']
#'             clonality = ifelse(cf.em >= (purity * 0.8), 'clonal', 'subclonal')
#'             pass = gene_level[which(gene_level$gene == gene), 'filter']
#'             
#'           } else {
#'             CnLR = NA
#'             cf.em = NA
#'             gene_cn = NA
#'             clonality = NA
#'             pass = NA
#'           }
#'           
#'         } else {
#'           CnLR = NA
#'           cf.em = NA
#'           gene_cn = NA
#'           clonality = NA
#'           pass = NA
#'         }
#'         
#'         
#'         #' t-test against baseline dipLogR
#'         # gene_chrom = as.integer(genes_hg19[which(genes_hg19$gene == gene), 'chrom'])
#'         # gene_start = as.integer(genes_hg19[which(genes_hg19$gene == gene), 'start'])
#'         # gene_end = as.integer(genes_hg19[which(genes_hg19$gene == gene), 'end'])
#'         # 
#'         # if(!is.na(CnLR) & CnLR > 0){
#'         #   one_sided_test = t.test(gene_snps, mu = fit$dipLogR, alternative = "greater")$p.value
#'         # } else if (!is.na(CnLR) & CnLR < 0) {
#'         #   one_sided_test = t.test(gene_snps, mu = fit$dipLogR, alternative = "less")$p.value
#'         # } else {
#'         #   one_sided_test = NA
#'         # }
#'         
#'         
#'         #' prepare output
#'         out = data.frame(gene = gene,
#'                          QC_all = facets_qc$facets_qc,
#'                          QC_minimum = qc_5parameters,
#'                          copy_state = gene_cn,
#'                          clonality = clonality,
#'                          cf.em = cf.em,
#'                          dipLogR_original = dipLogR_original,
#'                          cnlr = CnLR,
#'                          pass = pass, 
#'                          group = tumor_sample)
#'         
#'         gene_level_out = rbind(gene_level_out, out)
#'         facets_plots[[tumor_sample]] = all_plots
#'         
#'         rm(gene_snps, gene_chrom, gene_start, gene_end, gene_cn, 
#'            pass, clonality, cf.em, CnLR)
#' 
#'     })
#'   }
#'   rm(facets_qc, qc_5parameters)
#'   gene_level_out = gene_level_out[!duplicated(gene_level_out), ]
#' }
#' 
#' 
#' 
#' ##-----------------
#' ## manual inspection and re-run
#' ##-----------------
#' manual = multi_readSnpMatrix(filename = countMatrix_path, tumor_sample = 3)
#' fit = facetsSuite::run_facets(read_counts = manual, 
#'                               cval = cval,
#'                               min_nhet = min_het,
#'                               seed = seed,
#'                               genome = 'hg19')
#' fit$dipLogR
#' i = facetsSuite::cnlr_plot(fit, return_object = T)
#' ii = facetsSuite::valor_plot(fit, return_object = T)
#' iii = facetsSuite::icn_plot(fit, return_object = T)
#' iv = facetsSuite::cf_plot(fit, return_object = T)
#' dev.off()
#' i / ii/ iii/ iv
#' j = facets_fit_qc(fit)
#' j
#' 
#' 
#' samples_dipLogR = c(0.07104743)
#' 
#' 
#' ##-----------------
#' ## SECOND Run with proper dipLogR
#' ##-----------------
#' ## Parameters: (exclusively purity runs); not interested in gene_level alterations
#' cval = 100
#' seed = 100
#' min_het = 15
#' genome = 'hg19'
#' 
#' parameter_table$dipLogR = samples_dipLogR
#' dev.off()
#' gene_level_out = data.frame()
#' facets_plots = list()
#' parameter_table$ID = NA
#' for(tumor_sample in 1:nrow(parameter_table)){
#'   parameter_table$ID[tumor_sample] = sample_match$sample[grepl(parameter_table$name[tumor_sample], sample_match$path)]
#'   countMatrix = multi_readSnpMatrix(filename = countMatrix_path, 
#'                                     tumor_sample = parameter_table$tumor_sample[tumor_sample])
#'   
#'   fit = facetsSuite::run_facets(read_counts = countMatrix,
#'                                 cval = cval,
#'                                 min_nhet = min_het,
#'                                 genome = genome,
#'                                 seed = seed,
#'                                 dipLogR = parameter_table$dipLogR[tumor_sample])
#'   
#'   print(paste0(parameter_table$name[tumor_sample], ': facets_qc ', facets_fit_qc(fit)$facets_qc))
#'   filters = c('waterfall_filter_pass', 'hyper_seg_filter_pass', 
#'               'em_cncf_icn_discord_filter_pass', 'contamination_filter_pass')
#'   facets_qc = facets_fit_qc(fit)
#'   
#'   if(all(unlist(facets_qc[filters]))){
#'     qc_5parameters = TRUE
#'     print(paste0('min qc :', qc_5parameters))
#'   } else {
#'     qc_5parameters = FALSE
#'     print(paste0('min qc :', qc_5parameters))
#'   }
#'   
#'   i = facetsSuite::cnlr_plot(fit, return_object = T)
#'   ii = facetsSuite::valor_plot(fit, return_object = T)
#'   iii = facetsSuite::icn_plot(fit, return_object = T)
#'   iv = facetsSuite::cf_plot(fit, return_object = T)
#'   print(i / ii/ iii/ iv)
#'   all_plots = i / ii / iii / iv
#'   
#'   #' gene level changes
#'   gene_level = facetsSuite::gene_level_changes(fit,
#'                                                genome = 'hg19')
#'   dipLogR_original = fit$dipLogR
#'   print(paste0('Original dipLogR: ', dipLogR_original))
#'   purity = fit$purity
#'   
#'   #' run on clinically significant genes
#'   for(gene in c('CDKN2A',
#'                 'CDK4',
#'                 'CDK6',
#'                 'PTEN',
#'                 'EGFR',
#'                 'PDGFRA',
#'                 'KIT',
#'                 'KDR',
#'                 'MET',
#'                 'MDM2',
#'                 'MDM4',
#'                 'RB1',
#'                 'NF1',
#'                 'TP53')){
#'     try({
#'       #' Facets CnLR estimates
#'       if(facets_qc$facets_qc & qc_5parameters){
#'         if(gene %in% gene_level$gene){
#'           CnLR = gene_level[which(gene_level$gene == gene), 'median_cnlr_seg']
#'           CnLR = ifelse(is.nan(CnLR) | is.na(CnLR) | length(CnLR) == 0, NA, CnLR - dipLogR_original)
#'           cf.em = gene_level[which(gene_level$gene == gene), 'cf.em']
#'           cf.em = ifelse(is.nan(cf.em) | is.na(cf.em) | length(cf.em) == 0, NA, cf.em)
#'           gene_cn = gene_level[which(gene_level$gene == gene), 'cn_state']
#'           clonality = ifelse(cf.em >= (purity * 0.8), 'clonal', 'subclonal')
#'           pass = gene_level[which(gene_level$gene == gene), 'filter']
#'           
#'         } else {
#'           CnLR = NA
#'           cf.em = NA
#'           gene_cn = NA
#'           clonality = NA
#'           pass = NA
#'         }
#'         
#'       } else if (facets_qc$facets_qc | qc_5parameters) {
#'         if(gene %in% gene_level$gene){
#'           gene_chrom = as.integer(genes_hg19[which(genes_hg19$gene == gene), 'chrom'])
#'           gene_start = as.integer(genes_hg19[which(genes_hg19$gene == gene), 'start'])
#'           gene_end = as.integer(genes_hg19[which(genes_hg19$gene == gene), 'end'])
#'           gene_snps = fit$snps[which(fit$snps$chrom == gene_chrom & 
#'                                        fit$snps$maploc >= gene_start & 
#'                                        fit$snps$maploc <= gene_end), 'cnlr']
#'           gene_snps = as.numeric(gene_snps)
#'           CnLR = ifelse(is.nan(gene_snps) | is.na(gene_snps) | length(gene_snps) == 0, NA, mean(gene_snps, na.rm = T))
#'           cf.em = gene_level[which(gene_level$gene == gene), 'cf.em']
#'           cf.em = ifelse(is.nan(cf.em) | is.na(cf.em) | length(cf.em) == 0, NA, cf.em)
#'           gene_cn = gene_level[which(gene_level$gene == gene), 'cn_state']
#'           clonality = ifelse(cf.em >= (purity * 0.8), 'clonal', 'subclonal')
#'           pass = gene_level[which(gene_level$gene == gene), 'filter']
#'           
#'         } else {
#'           CnLR = NA
#'           cf.em = NA
#'           gene_cn = NA
#'           clonality = NA
#'           pass = NA
#'         }
#'         
#'       } else {
#'         CnLR = NA
#'         cf.em = NA
#'         gene_cn = NA
#'         clonality = NA
#'         pass = NA
#'       }
#'       
#'       
#'       #' t-test against baseline dipLogR
#'       # gene_chrom = as.integer(genes_hg19[which(genes_hg19$gene == gene), 'chrom'])
#'       # gene_start = as.integer(genes_hg19[which(genes_hg19$gene == gene), 'start'])
#'       # gene_end = as.integer(genes_hg19[which(genes_hg19$gene == gene), 'end'])
#'       # 
#'       # if(!is.na(CnLR) & CnLR > 0){
#'       #   one_sided_test = t.test(gene_snps, mu = fit$dipLogR, alternative = "greater")$p.value
#'       # } else if (!is.na(CnLR) & CnLR < 0) {
#'       #   one_sided_test = t.test(gene_snps, mu = fit$dipLogR, alternative = "less")$p.value
#'       # } else {
#'       #   one_sided_test = NA
#'       # }
#'       
#'       
#'       #' prepare output
#'       out = data.frame(gene = gene,
#'                        QC_all = facets_qc$facets_qc,
#'                        QC_minimum = qc_5parameters,
#'                        copy_state = gene_cn,
#'                        clonality = clonality,
#'                        cf.em = cf.em,
#'                        dipLogR_original = dipLogR_original,
#'                        cnlr = CnLR,
#'                        pass = pass, 
#'                        group = tumor_sample)
#'       
#'       gene_level_out = rbind(gene_level_out, out)
#'       facets_plots[[tumor_sample]] = all_plots
#'       
#'       rm(gene_snps, gene_chrom, gene_start, gene_end, gene_cn, 
#'          pass, clonality, cf.em, CnLR)
#'       
#'     })
#'   }
#'   
#'   rm(facets_qc, qc_5parameters)
#'   create_facets_output(facets_output = fit, 
#'                        directory = paste0(getwd(), '/', gsub("/.*$", "", countMatrix_path), '/'), 
#'                        sample_id = parameter_table$ID[tumor_sample])
#' }
#' 
#' 
#' 
#' 
#' rm(countMatrix, countMatrix_raw, gene_level, ii, iii, iv, out, i, filters, 
#'    tumor_sample, seed, samples_dipLogR, samples, purity, gene, cval, all_plots, facets_plots, fit)
#' 
#' gene_level_out$name = NA
#' for(i in unique(gene_level_out$group)){
#'   gene_level_out$name[which(gene_level_out$group == i)] = parameter_table$name[i]
#' }
#' 
#' write.table(x = gene_level_out, file = paste0(gsub("/.*$", "", countMatrix_path), '/', 'gene_level_out.txt'), sep = '\t', quote = F, row.names = F)
#' write.table(parameter_table, file = paste0(gsub("/.*$", "", countMatrix_path), '/', 'parameterTable.txt'), sep = '\t', quote = F, row.names = F)
#' 
#' 
#' 
#' ##' test:
#' rds = load('C-000499/s_C_000499_L002_d/s_C_000499_L002_d.Rdata')  

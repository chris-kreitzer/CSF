## Automated pipeline for copy number alteration detection
## using Facets!
## 
## start: 09/11/2022
## chris-kreitzer

clean()
gc()
.rs.restartR()
setup(working.path = '~/Documents/MSKCC/Subhi/CSF/')
library(patchwork)
source('~/Documents/GitHub/CSF/Scripts/UtilityFunctions.R')

## C_006884: countMatrix
countMatrix_path = 'C_THJ/C_THJ_countMatrix.dat.gz'
countMatrix_raw = read.csv(file = countMatrix_path, sep = ',')
samples = grep(pattern = 'File*', colnames(countMatrix_raw))
samples = (length(samples) - 4) / 4

## Parameters: (exclusively purity runs); not interested in gene_level alterations
cval = 100
seed = 100
min_het = 15
genome = 'hg19'

##-----------------
## First RUN
##-----------------
parameter_table = data.frame(tumor_sample = c(2,3,4,5),
                             name = c('P-0009511-T02-IM5',
                                      'P-0009511-T03-IM6',
                                      's_C_THJT7J_S001_d',
                                      's_C_THJT7J_S002_d'),
                             dipLogR = NA)
gene_level_out = data.frame()
facets_plots = list()
for(tumor_sample in 1:nrow(parameter_table)){
  countMatrix = multi_readSnpMatrix(filename = countMatrix_path, 
                                    tumor_sample = parameter_table$tumor_sample[tumor_sample])
  
  fit = facetsSuite::run_facets(read_counts = countMatrix,
                                cval = cval,
                                min_nhet = min_het,
                                genome = genome,
                                seed = seed)
  
  print(paste0(parameter_table$name[tumor_sample], ': facets_qc ', facets_fit_qc(fit)$facets_qc))
  filters = c('waterfall_filter_pass', 'hyper_seg_filter_pass', 
              'em_cncf_icn_discord_filter_pass', 'contamination_filter_pass')
  facets_qc = facets_fit_qc(fit)
  
  if(all(unlist(facets_qc[filters]))){
    qc_5parameters = TRUE
    print(paste0('min qc :', qc_5parameters))
  } else {
    qc_5parameters = FALSE
    print(paste0('min qc :', qc_5parameters))
  }
  
  i = facetsSuite::cnlr_plot(fit, return_object = T)
  ii = facetsSuite::valor_plot(fit, return_object = T)
  iii = facetsSuite::icn_plot(fit, return_object = T)
  iv = facetsSuite::cf_plot(fit, return_object = T)
  print(i / ii/ iii/ iv)
  all_plots = i / ii / iii / iv
  
  #' gene level changes
  gene_level = facetsSuite::gene_level_changes(fit,
                                               genome = 'hg19')
  dipLogR_original = fit$dipLogR
  print(paste0('Original dipLogR: ', dipLogR_original))
  purity = fit$purity
  
  #' run on clinically significant genes
  for(gene in c('CDKN2A',
                'CDK4',
                'CDK6',
                'PTEN',
                'EGFR',
                'PDGFRA',
                'KIT',
                'KDR',
                'MET',
                'MDM2',
                'MDM4',
                'RB1',
                'NF1',
                'TP53')){
    
    #' Facets CnLR estimates
    if(facets_qc$facets_qc & qc_5parameters){
     CnLR = gene_level[which(gene_level$gene == gene), 'median_cnlr_seg']
     CnLR = CnLR - dipLogR_original
     cf.em = gene_level[which(gene_level$gene == gene), 'cf.em']
     
    } else if (facets_qc$facets_qc | qc_5parameters) {
      gene_chrom = as.integer(genes_hg19[which(genes_hg19$gene == gene), 'chrom'])
      gene_start = as.integer(genes_hg19[which(genes_hg19$gene == gene), 'start'])
      gene_end = as.integer(genes_hg19[which(genes_hg19$gene == gene), 'end'])
      gene_snps = fit$snps[which(fit$snps$chrom == gene_chrom & 
                                   fit$snps$maploc >= gene_start & 
                                   fit$snps$maploc <= gene_end), 'cnlr']
      gene_snps = as.numeric(gene_snps)
      CnLR = mean(gene_snps)
      cf.em = gene_level[which(gene_level$gene == gene), 'cf.em']
    } else {
      CnLR = NA
      cf.em = NA
    }
    
    #' assign clonality
    clonality = ifelse(cf.em >= (purity * 0.8), 'clonal', 'subclonal')
    gene_cn = gene_level[which(gene_level$gene == gene), 'cn_state']
    pass = gene_level[which(gene_level$gene == gene), 'filter']
    
    #' t-test against baseline dipLogR
    gene_chrom = as.integer(genes_hg19[which(genes_hg19$gene == gene), 'chrom'])
    gene_start = as.integer(genes_hg19[which(genes_hg19$gene == gene), 'start'])
    gene_end = as.integer(genes_hg19[which(genes_hg19$gene == gene), 'end'])
    
    #' snps
    gene_snps = fit$snps[which(fit$snps$chrom == gene_chrom & 
                                 fit$snps$maploc >= gene_start & 
                                 fit$snps$maploc <= gene_end), 'cnlr']
    gene_snps = as.numeric(gene_snps)
    
    if(CnLR > 0){
      one_sided_test = t.test(gene_snps, mu = fit$dipLogR, alternative = "greater")$p.value
    } else {
      one_sided_test = t.test(gene_snps, mu = fit$dipLogR, alternative = "less")$p.value
    }
    
    
    #' prepare output
    out = data.frame(gene = gene,
                     QC_all = facets_qc$facets_qc,
                     QC_minimum = qc_5parameters,
                     copy_state = gene_cn,
                     clonality = clonality,
                     cf.em = cf.em,
                     dipLogR_original = dipLogR_original,
                     cnlr = CnLR,
                     one_sided_pvalue = one_sided_test,
                     pass = pass, 
                     group = tumor_sample)
    
    gene_level_out = rbind(gene_level_out, out)
    facets_plots[[tumor_sample]] = all_plots
    
    rm(one_sided_test, gene_snps, gene_chrom, gene_start, gene_end, gene_cn, 
       pass, clonality, cf.em, CnLR)
  }
  rm(facets_qc, qc_5parameters)
}


##-----------------
## manual inspection and re-run
##-----------------
manual = multi_readSnpMatrix(filename = 'C_006884/C_006884_countmatrix_dat.gz', tumor_sample = 6)
fit = facetsSuite::run_facets(read_counts = manual, 
                              cval = cval,
                              min_nhet = min_het,
                              seed = seed,
                              genome = 'hg19')
fit$dipLogR
i = facetsSuite::cnlr_plot(fit, return_object = T)
ii = facetsSuite::valor_plot(fit, return_object = T)
iii = facetsSuite::icn_plot(fit, return_object = T)
iv = facetsSuite::cf_plot(fit, return_object = T)
i / ii/ iii/ iv

j = facets_fit_qc(fit)
j
fit$segs

samples_dipLogR = c(-0.09901, 0.0934804, -0.02306294663361, 0.012657054)

  
##-----------------
## SECOND Run with proper dipLogR
##-----------------
## Parameters: (exclusively purity runs); not interested in gene_level alterations
cval = 100
seed = 100
min_het = 15
genome = 'hg19'

parameter_table$dipLogR = samples_dipLogR

gene_level_out = data.frame()
facets_plots = list()
for(tumor_sample in 1:nrow(parameter_table)){
  countMatrix = multi_readSnpMatrix(filename = countMatrix_path, 
                                    tumor_sample = parameter_table$tumor_sample[tumor_sample])
  
  fit = facetsSuite::run_facets(read_counts = countMatrix,
                                cval = cval,
                                min_nhet = min_het,
                                genome = genome,
                                seed = seed,
                                dipLogR = parameter_table$dipLogR[tumor_sample])
  
  print(paste0(parameter_table$name[tumor_sample], ': facets_qc ', facets_fit_qc(fit)$facets_qc))
  filters = c('waterfall_filter_pass', 'hyper_seg_filter_pass', 
              'em_cncf_icn_discord_filter_pass', 'contamination_filter_pass')
  facets_qc = facets_fit_qc(fit)
  
  if(all(unlist(facets_qc[filters]))){
    qc_5parameters = TRUE
    print(paste0('min qc :', qc_5parameters))
  } else {
    qc_5parameters = FALSE
    print(paste0('min qc :', qc_5parameters))
  }
  
  i = facetsSuite::cnlr_plot(fit, return_object = T)
  ii = facetsSuite::valor_plot(fit, return_object = T)
  iii = facetsSuite::icn_plot(fit, return_object = T)
  iv = facetsSuite::cf_plot(fit, return_object = T)
  print(i / ii/ iii/ iv)
  all_plots = i / ii / iii / iv
  
  #' gene level changes
  gene_level = facetsSuite::gene_level_changes(fit,
                                               genome = 'hg19')
  dipLogR_original = fit$dipLogR
  print(paste0('Original dipLogR: ', dipLogR_original))
  purity = fit$purity
  
  #' run on clinically significant genes
  for(gene in c('CDKN2A',
                'CDK4',
                'CDK6',
                'PTEN',
                'EGFR',
                'PDGFRA',
                'KIT',
                'KDR',
                'MET',
                'MDM2',
                'MDM4',
                'RB1',
                'NF1',
                'TP53')){
    
    #' Facets CnLR estimates
    if(facets_qc$facets_qc & qc_5parameters){
      CnLR = gene_level[which(gene_level$gene == gene), 'median_cnlr_seg']
      CnLR = CnLR - dipLogR_original
      cf.em = gene_level[which(gene_level$gene == gene), 'cf.em']
      
    } else if (facets_qc$facets_qc | qc_5parameters) {
      gene_chrom = as.integer(genes_hg19[which(genes_hg19$gene == gene), 'chrom'])
      gene_start = as.integer(genes_hg19[which(genes_hg19$gene == gene), 'start'])
      gene_end = as.integer(genes_hg19[which(genes_hg19$gene == gene), 'end'])
      gene_snps = fit$snps[which(fit$snps$chrom == gene_chrom & 
                                   fit$snps$maploc >= gene_start & 
                                   fit$snps$maploc <= gene_end), 'cnlr']
      gene_snps = as.numeric(gene_snps)
      CnLR = mean(gene_snps)
      cf.em = gene_level[which(gene_level$gene == gene), 'cf.em']
    } else {
      CnLR = NA
      cf.em = NA
    }
    
    #' assign clonality
    clonality = ifelse(cf.em >= (purity * 0.8), 'clonal', 'subclonal')
    gene_cn = gene_level[which(gene_level$gene == gene), 'cn_state']
    pass = gene_level[which(gene_level$gene == gene), 'filter']
    
    #' t-test against baseline dipLogR
    gene_chrom = as.integer(genes_hg19[which(genes_hg19$gene == gene), 'chrom'])
    gene_start = as.integer(genes_hg19[which(genes_hg19$gene == gene), 'start'])
    gene_end = as.integer(genes_hg19[which(genes_hg19$gene == gene), 'end'])
    
    #' snps
    gene_snps = fit$snps[which(fit$snps$chrom == gene_chrom & 
                                 fit$snps$maploc >= gene_start & 
                                 fit$snps$maploc <= gene_end), 'cnlr']
    gene_snps = as.numeric(gene_snps)
    
    if(CnLR > 0){
      one_sided_test = t.test(gene_snps, mu = fit$dipLogR, alternative = "greater")$p.value
    } else {
      one_sided_test = t.test(gene_snps, mu = fit$dipLogR, alternative = "less")$p.value
    }
    
    
    #' prepare output
    out = data.frame(gene = gene,
                     QC_all = facets_qc$facets_qc,
                     QC_minimum = qc_5parameters,
                     copy_state = gene_cn,
                     clonality = clonality,
                     cf.em = cf.em,
                     dipLogR_original = dipLogR_original,
                     cnlr = CnLR,
                     one_sided_pvalue = one_sided_test,
                     pass = pass, 
                     group = tumor_sample)
    
    gene_level_out = rbind(gene_level_out, out)
    facets_plots[[tumor_sample]] = all_plots
    
    rm(one_sided_test, gene_snps, gene_chrom, gene_start, gene_end, gene_cn, 
       pass, clonality, cf.em, CnLR)
  }
  rm(facets_qc, qc_5parameters)
  create_facets_output(facets_output = fit, directory = paste0(getwd(), '/C_THJ/'), sample_id = parameter_table$name[tumor_sample])
}



gene_level_out$name = NA
for(i in unique(gene_level_out$group)){
  gene_level_out$name[which(gene_level_out$group == i)] = parameter_table$name[i]
}

write.table(x = gene_level_out, file = 'C_THJ/gene_level_out.txt', sep = '\t', quote = F, row.names = F)
  
  
  
  
  
  












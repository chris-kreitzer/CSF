##-----------------
## C_006884:

clean()
gc()
source('~/Documents/GitHub/CSF/Scripts/UtilityFunctions.R')
library(patchwork)

gene_level_out = data.frame()
for(tumor_sample in c(2, 3, 4, 5, 6)){
  #rm(fit, gene_level, gene_level_out)
  countMatrix = multi_readSnpMatrix(filename = 'C_006884/C_006884_countmatrix_dat.gz', tumor_sample = tumor_sample)
  
  fit = facetsSuite::run_facets(read_counts = countMatrix,
                                cval = 100, 
                                genome = 'hg19',
                                seed = 100)
                                
  print(paste0(tumor_sample, ': facets_qc ', facets_fit_qc(fit)$facets_qc))
  
  i = facetsSuite::cnlr_plot(fit, return_object = T)
  ii = facetsSuite::valor_plot(fit, return_object = T)
  iii = facetsSuite::icn_plot(fit, return_object = T)
  iv = facetsSuite::cf_plot(fit, return_object = T)
  print(i / ii/ iii/ iv)
  
  #' gene level
  gene_level = facetsSuite::gene_level_changes(fit,
                                               genome = 'hg19')
  dipLogR_original = fit$dipLogR
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
    
    median_cnlr_observed = gene_level[which(gene_level$gene == gene), 'median_cnlr_seg']
    median_cnlr_observed = median_cnlr_observed - dipLogR_original
    cf.em = gene_level[which(gene_level$gene == gene), 'cf.em']
    
    #' assign clonality
    clonality = ifelse(cf.em >= (purity * 0.8), 'clonal', 'subclonal')
    gene_cn = gene_level[which(gene_level$gene == gene), 'cn_state']
    pass = gene_level[which(gene_level$gene == gene), 'filter']
    
    out = data.frame(gene = gene,
                     cnlr = median_cnlr_observed,
                     cf.em = cf.em,
                     dipLogR_original = dipLogR_original,
                     copy_state = gene_cn,
                     clonality = clonality,
                     pass = pass, 
                     group = tumor_sample)
    
    gene_level_out = rbind(gene_level_out, out)
  }
}

#' add cBIO annotation
cBIO = gene_level_out[1:14, ]
cBIO$group = 'cBIO-DMP1'
cBIO$copy_state[which(cBIO$gene == 'CDK6')] = 'DIPLOID'
cBIO$copy_state[which(cBIO$gene == 'EGFR')] = 'DIPLOID'
cBIO$copy_state[which(cBIO$gene == 'MET')] = 'DIPLOID'
cBIO$clonality[which(cBIO$gene == 'CDK6')] = 'clonal'
cBIO2 = cBIO
cBIO2$group = 'cBIO-DMP2'

gene_level_out = rbind(cBIO, cBIO2, gene_level_out)
gene_level_out$group[which(gene_level_out$group == 2)] = 'Facets-DMP1'
gene_level_out$group[which(gene_level_out$group == 3)] = 'Facets-DMP2'
gene_level_out$group[which(gene_level_out$group == 4)] = 'CSF1'
gene_level_out$group[which(gene_level_out$group == 5)] = 'CSF2'
gene_level_out$group[which(gene_level_out$group == 6)] = 'CSF3'
gene_level_out$group = factor(gene_level_out$group, levels = c('cBIO-DMP1', 'Facets-DMP1', 'cBIO-DMP2', 'Facets-DMP2', 'CSF1', 'CSF2', 'CSF3'))
gene_level_out$clonality = ifelse(gene_level_out$clonality == 'clonal', '', '*')
gene_level_out$copy_state[which(gene_level_out$copy_state == 'CNLOH')] = 'LOSS'
gene_level_out$copy_state[which(gene_level_out$copy_state == 'HOMDEL')] = 'LOSS'
gene_level_out$copy_state[which(gene_level_out$copy_state == 'HETLOSS')] = 'LOSS'
gene_level_out$copy_state[which(gene_level_out$copy_state == 'GAIN (many states)')] = 'GAIN'

DMP_C006884 = ggplot(gene_level_out, aes(x = group, y = gene, fill = copy_state)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) +
  scale_fill_manual(values = c('LOSS' = 'lightblue',
                               'DIPLOID' = 'grey85',
                               'GAIN' = 'firebrick3'),
                    name = 'CNA') +
  coord_fixed() +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  geom_text(aes(label = clonality), color = "white", size = 4) +
  theme_std(base_size = 14) +
  theme(axis.ticks = element_blank(), 
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(angle = 65, vjust = 0, hjust = 0))

DMP_C006884
ggsave(filename = 'C_006884/geneLevel_DMPnormal.pdf', plot = DMP_C006884, device = 'pdf', width = 4, height = 8)


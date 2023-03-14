##----------------+
## Cohort overview CSF
## General statistics
## about CSF and Tumors
##----------------+

## start: 03/12/2023
## revision: 03/13/2023
## 
## chris-kreitzer


clean()
setwd('~/Documents/MSKCC/11_CSF/')
source('~/Documents/GitHub/DryClean_Facets/Scripts/CnLR_plot.R')
source('~/Documents/GitHub/Y_chromosome_loss/PanCancer/Scripts/UtilityFunctions.R')
source('~/Documents/MSKCC/11_CSF/02_Scripts/CnLR_plot.R')

library(cowplot)
library(patchwork)
library(ggpubr)

database = readxl::read_excel('00_Data/Database_Final_march12.xlsx', sheet = 1)
#database = database[which(database$Fit %in% c('full', 'reduced')), ]
#colnames(database)[1] = 'Sample.ID'
all_files = list.files('01_countmatrices/', full.names = T, recursive = T)
copynumberstates = facetsSuite:::copy_number_states
GOIs = c('CDKN2A','CDK4','CDK6','PTEN','EGFR','RET',
         'GATA3','PDGFRA','KIT','KDR','MET','MDM2',
         'MDM4','RB1','TP53')
pairing = read.csv('PairedAnalysis_031023.txt', sep = '\t')
pairing = merge(pairing, database[,c('Sample.ID', 'Fit', 'TYPE')], by = 'Sample.ID', all.x = T)


##----------------+
## Mutations overview
##----------------+
data_mutations = rbind(database[which(database$TYPE == 'TUMOR'), ],
                       database[which(database$TYPE == 'CSF' & database$CSF_PASS == 'Pass'), ])

mutations_out = data.frame()
for(i in 1:nrow(data_mutations)){
  print(i)
  id = data_mutations$Sample.ID[i]
  oncogenic = data_mutations$Oncogenic_Muts_2023[i]
  oncogenic = ifelse(is.na(oncogenic), 0, 
                     length(unlist(strsplit(x = oncogenic, split = ','))))
  
  non_oncogenic = data_mutations$Other_non_oncogenic_Muts[i]
  non_oncogenic = ifelse(is.na(non_oncogenic), 0, 
                         length(unlist(strsplit(x = non_oncogenic, split = ','))))
  
  type = data_mutations$TYPE[i]
  
  out = data_frame(id = id,
                   type = type,
                   oncogenic = oncogenic,
                   non_oncogenic = non_oncogenic)
  mutations_out = rbind(mutations_out, out)
}

mutations_out$type = factor(mutations_out$type, levels = c('TUMOR', 'CSF'))


##-- Visualization
nMutations = ggplot(mutations_out, aes(x = type, y = oncogenic)) +
  geom_boxplot(outlier.colour = 'white') +
  geom_jitter(width = 0.2, size = 0.5) +
  theme_std(base_size = 14, base_line_size = 0.1) +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0, 30)) +
  panel_border(size = 2, color = 'black') +
  theme(aspect.ratio = 1.5,
        axis.ticks.y = element_line(linewidth = 0.85)) +
  labs(y = '# Mutations oncogenic', x = '') 

nMutations = nMutations + annotate(geom = 'text', x = 1.5, y = 28, label = 'P=0.00143')


#' non-oncogenic mutations
nononc_Mutations = ggplot(mutations_out, aes(x = type, y = log10(non_oncogenic + 1))) +
  geom_boxplot(outlier.colour = 'white') +
  geom_jitter(width = 0.2, size = 0.5) +
  theme_std(base_size = 14, base_line_size = 0.1) +
  scale_y_continuous(expand = c(0,0)) +
  panel_border(size = 2, color = 'black') +
  theme(aspect.ratio = 1.5,
        axis.ticks.y = element_line(linewidth = 0.85)) +
  labs(y = '# Mutations all [log10]', x = '')


nononc_Mutations = nononc_Mutations + annotate(geom = 'text', x = 1.5, y = 2.5, label = 'P=0.117')

mutations = nMutations + nononc_Mutations
ggsave_golden(filename = '05_Plots/Mutations_TumorCSF.pdf', plot = mutations, width = 8)




##----------------+
## Alteration OVERVIEW ALL
##----------------+
all_out = data.frame()
for(i in 1:nrow(database)){
  try({
    print(i)
    id = database$Sample.ID[i]
    fit = database$Fit[i]
    
    if(fit == 'full'){
      file = all_files[grep(pattern = id, x = all_files)]
      filerds = file[grep(pattern = '.rds$', x = file)]
      data = readRDS(filerds)
      
      ##-- run Tests
      qc = facets_fit_qc(facets_output = data)
      gene = facetsSuite::gene_level_changes(facets_output = data, genome = 'hg19', algorithm = 'em')
      gene$string = paste(qc$wgd, gene$tcn.em, gene$mcn, gene$lcn.em, sep = ':')
      gene = merge(gene, copynumberstates[,c('map_string', 'numeric_call')],  by.x = 'string', by.y = 'map_string', all.x = T)
      
      fga = qc$fga
      purity = qc$purity
      diplogr = qc$dipLogR
      n_segs = qc$n_segs
      frac_below_diplogr = qc$frac_below_dipLogR
      waterfall = qc$sd_cnlr_residual
      qc_facets = qc$facets_qc
      
      ##-- alterations
      table_out = table(gene$numeric_call)
      
      frac_loss = table_out[which(names(table_out) %in% c('-1', '-2'))]
      frac_loss = ifelse(length(frac_loss) == 2, (frac_loss[[1]] + frac_loss[[2]]), frac_loss[[1]])
      
      frac_gain = table_out[which(names(table_out) %in% c('1', '2'))]
      frac_gain = ifelse(length(frac_gain) == 2, (frac_gain[[1]] + frac_gain[[2]]), frac_gain[[1]])
      
      out = data.frame(id = id,
                       fga = fga,
                       purity = purity,
                       diplogr = diplogr,
                       n_segs = n_segs,
                       frac_below_diplogr = frac_below_diplogr,
                       waterfall_sd = waterfall,
                       qc = qc_facets,
                       frac_loss = (frac_loss) / sum(table(gene$numeric_call)),
                       frac_gain = (frac_gain) / sum(table(gene$numeric_call)))
      
    } else if (fit == 'reduced'){
      file = all_files[grep(pattern = id, x = all_files)]
      file = file[grep(pattern = 'adjusted.*', x = file)]
      filerds = file[grep(pattern = 'rds$', file)]
      data = readRDS(filerds)
      
      ##-- run Tests
      qc = facets_fit_qc(facets_output = data)
      gene = facetsSuite::gene_level_changes(facets_output = data, genome = 'hg19', algorithm = 'em')
      gene$string = paste(qc$wgd, gene$tcn.em, gene$mcn, gene$lcn.em, sep = ':')
      gene = merge(gene, copynumberstates[,c('map_string', 'numeric_call')],  by.x = 'string', by.y = 'map_string', all.x = T)
      
      fga = qc$fga
      purity = qc$purity
      diplogr = qc$dipLogR
      n_segs = qc$n_segs
      frac_below_diplogr = qc$frac_below_dipLogR
      waterfall = qc$sd_cnlr_residual
      qc_facets = qc$facets_qc
      
      ##-- alterations
      table_out = table(gene$numeric_call)
      
      frac_loss = table_out[which(names(table_out) %in% c('-1', '-2'))]
      frac_loss = ifelse(length(frac_loss) == 2, (frac_loss[[1]] + frac_loss[[2]]), frac_loss[[1]])
      
      frac_gain = table_out[which(names(table_out) %in% c('1', '2'))]
      frac_gain = ifelse(length(frac_gain) == 2, (frac_gain[[1]] + frac_gain[[2]]), frac_gain[[1]])
      
      out = data.frame(id = id,
                       fga = fga,
                       purity = purity,
                       diplogr = diplogr,
                       n_segs = n_segs,
                       frac_below_diplogr = frac_below_diplogr,
                       waterfall_sd = waterfall,
                       qc = qc_facets,
                       frac_loss = (frac_loss) / sum(table(gene$numeric_call)),
                       frac_gain = (frac_gain) / sum(table(gene$numeric_call)))
      
    } else next
    
    all_out = rbind(all_out, out)
  })
}


all_out = merge(all_out, database_short, by.x = 'id', by.y = 'Sample.ID', all.x = T)
all_out = merge(all_out, database[,c('Sample.ID', 'TYPE')], by.x = 'id', by.y = 'Sample.ID', all.x = T)




##-- make some plots
full_fit = all_out[which(all_out$Fit == 'full'), ]
full_fit$TYPE = factor(full_fit$TYPE, levels = c('TUMOR', 'CSF'))


FGA = ggplot(full_fit, aes(x = TYPE, y = fga)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 0.5) +
  theme_std(base_size = 14, base_line_size = 0.1) +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0, 1)) +
  panel_border(size = 2, color = 'black') +
  theme(aspect.ratio = 1.5,
        axis.ticks.y = element_line(linewidth = 0.85)) +
  labs(y = 'Fraction Genome Altered', x = '')

ggsave_golden(filename = '05_Plots/FGA_TumorCSF.pdf', plot = FGA, width = 6)


##-- Purity
Purity = ggplot(full_fit, aes(x = TYPE, y = purity)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 0.5) +
  theme_std(base_size = 14, base_line_size = 0.1) +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0, 1),
                     labels = c('0%', '25%', '50%', '75%', '100%')) +
  panel_border(size = 2, color = 'black') +
  theme(aspect.ratio = 1.5,
        axis.ticks.y = element_line(linewidth = 0.85)) +
  labs(y = 'Purity', x = '')

ggsave_golden(filename = '05_Plots/Purity_TumorCSF.pdf', plot = Purity, width = 6)


frac_gain_loss = full_fit[,c('TYPE', 'frac_loss', 'frac_gain')]
a = data.frame(type = frac_gain_loss$TYPE, 
               what = 'loss',
               value = frac_gain_loss$frac_loss)

b = data.frame(type = frac_gain_loss$TYPE, 
               what = 'gain',
               value = frac_gain_loss$frac_gain)

frac_gain_loss = rbind(a,b)

frac_gain_loss_plot = ggplot(frac_gain_loss, aes(x = type, y = value, color = what)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.45) +
  geom_jitter(position = position_dodge(width = 0.8), size = 0.45) +
  theme_std(base_size = 14, base_line_size = 0.1) +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0, 1),
                     labels = c('0%', '25%', '50%', '75%', '100%')) +
  scale_color_manual(values = c('gain' = 'red3',
                                'loss' = 'blue4'), name = '') +
  panel_border(size = 2, color = 'black') +
  theme(aspect.ratio = 1.5,
        axis.ticks.y = element_line(linewidth = 0.85),
        legend.position = 'top') +
  labs(y = 'Percent Gain and Loss', x = '')

ggsave_golden(filename = '05_Plots/Fraction_Gain_Loss.pdf', plot = frac_gain_loss_plot, width = 6)



##-- frac_below_diplogr
fraction_below_diplogr = ggplot(full_fit, aes(x = TYPE, y = frac_below_diplogr)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 0.5) +
  theme_std(base_size = 14, base_line_size = 0.1) +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0, 1)) +
  panel_border(size = 2, color = 'black') +
  theme(aspect.ratio = 1.5,
        axis.ticks.y = element_line(linewidth = 0.85)) +
  labs(y = 'Fraction below dipLogR', x = '')

ggsave_golden(filename = '05_Plots/Fraction_below_diplogr.pdf', plot = fraction_below_diplogr, width = 6)


##-- waterfall plot
waterfall_plot = ggplot(full_fit, aes(x = TYPE, y = waterfall_sd)) +
  geom_boxplot(outlier.size = 0.5) +
  geom_jitter(width = 0.2, size = 0.5) +
  theme_std(base_size = 14, base_line_size = 0.1) +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0, 1.5)) +
  panel_border(size = 2, color = 'black') +
  theme(aspect.ratio = 1.5,
        axis.ticks.y = element_line(linewidth = 0.85)) +
  labs(y = 'Log2 Ratio residual [sd]', x = '')


ggsave_golden(filename = '05_Plots/Waterfall.pdf', plot = waterfall_plot, width = 6)






##----------------+
## CDKN2A; paired analysis
##----------------+
pairs = read.csv('PairedAnalysis_031023.txt', sep = '\t')
pairs_use = database[which(database$Sample.ID %in% pairs$Sample.ID[which(pairs$Paired_Analysis == 'use')]), c('Sample.ID', 'Patient.ID', 'TYPE', 'Fit') ]

cdk_pairs = cd[which(cd$Patient.ID %in% c('C-000499',
                                          'C-06J46F',
                                          'C-09T7CL',
                                          'C-56D79J',
                                          'C-A32F9N')), ]

paired_analysis = data.frame()
for(i in 1:nrow(pairs_use)){
  try({
    print(i)
    id = pairs_use$Sample.ID[i]
    p.id = unique(pairs_use$Patient.ID[i])
    fit = pairs_use$Fit[i]
    
    if(fit == 'full'){
      file = all_files[grep(pattern = id, x = all_files)]
      filerds = file[grep(pattern = '.rds$', x = file)]
      data = readRDS(filerds)
      
      ##-- run Tests
      qc = facets_fit_qc(facets_output = data)
      gene = facetsSuite::gene_level_changes(facets_output = data, genome = 'hg19', algorithm = 'em')
      gene$string = paste(qc$wgd, gene$tcn.em, gene$mcn, gene$lcn.em, sep = ':')
      gene = merge(gene, copynumberstates[,c('map_string', 'numeric_call')],  by.x = 'string', by.y = 'map_string', all.x = T)
      cdkn2a = gene[which(gene$gene == 'CDKN2A'), c('gene', 'seg_start', 'seg_end', 'cn_state', 'numeric_call')]
      cdkn2a$id = id
      cdkn2a$p.id = p.id
      
      fga = qc$fga
      purity = qc$purity
      diplogr = qc$dipLogR
      n_segs = qc$n_segs
      frac_below_diplogr = qc$frac_below_dipLogR
      waterfall = qc$sd_cnlr_residual
      qc_facets = qc$facets_qc
      
    } else {
      file = all_files[grep(pattern = id, x = all_files)]
      file = file[grep(pattern = 'adjusted.*', x = file)]
      filerds = file[grep(pattern = 'rds$', file)]
      data = readRDS(filerds)
      
      ##-- run Tests
      qc = facets_fit_qc(facets_output = data)
      gene = facetsSuite::gene_level_changes(facets_output = data, genome = 'hg19', algorithm = 'em')
      gene$string = paste(qc$wgd, gene$tcn.em, gene$mcn, gene$lcn.em, sep = ':')
      gene = merge(gene, copynumberstates[,c('map_string', 'numeric_call')],  by.x = 'string', by.y = 'map_string', all.x = T)
      cdkn2a = gene[which(gene$gene == 'CDKN2A'), c('gene', 'seg_start', 'seg_end', 'cn_state', 'numeric_call')]
      cdkn2a$id = id
      cdkn2a$p.id = p.id
      
      fga = qc$fga
      purity = qc$purity
      diplogr = qc$dipLogR
      n_segs = qc$n_segs
      frac_below_diplogr = qc$frac_below_dipLogR
      waterfall = qc$sd_cnlr_residual
      qc_facets = qc$facets_qc
    }
  })
  paired_analysis = rbind(paired_analysis, cdkn2a)
}
    


##----------------+
## SAME CALLS FOR TUMOR AND CSF      
##----------------+
same_cdkn2a = c('C-001521', 'C-002114', 'C-06J46F', 'C-A32F9N', 'C-D2WW09', 'C-NRMAAD', 'C-RT2H2R')
cdkn2a_pairs = pairs[which(pairs$Patient.ID %in% same_cdkn2a & pairs$Paired_Analysis == 'use'), ]
cdkn2a_pairs = merge(cdkn2a_pairs, database[,c('Sample.ID', 'Fit', 'TYPE')], by = 'Sample.ID', all.x = T)

cdkn2a = data.frame()
for(i in 1:nrow(cdkn2a_pairs)){
  try({
    print(i)
    id = cdkn2a_pairs$Sample.ID[i]
    fit = cdkn2a_pairs$Fit[i]
    p.id = unique(cdkn2a_pairs$Patient.ID[i])
    type = unique(cdkn2a_pairs$TYPE[i])
    
    if(fit == 'full'){
      file = all_files[grep(pattern = id, x = all_files)]
      filerds = file[grep(pattern = '.rds$', x = file)]
      data = readRDS(filerds)
      
      ##-- run Tests
      qc = facets_fit_qc(facets_output = data)
      gene = facetsSuite::gene_level_changes(facets_output = data, genome = 'hg19', algorithm = 'em')
      gene$string = paste(qc$wgd, gene$tcn.em, gene$mcn, gene$lcn.em, sep = ':')
      gene = merge(gene, copynumberstates[,c('map_string', 'numeric_call')],  by.x = 'string', by.y = 'map_string', all.x = T)
      cdkn2a_call = gene$numeric_call[which(gene$gene == 'CDKN2A')]
      
      fga = qc$fga
      purity = qc$purity
      diplogr = qc$dipLogR
      n_segs = qc$n_segs
      frac_below_diplogr = qc$frac_below_dipLogR
      waterfall = qc$sd_cnlr_residual
      qc_facets = qc$facets_qc
      
      ##-- alterations
      table_out = table(gene$numeric_call)
      
      frac_loss = table_out[which(names(table_out) %in% c('-1', '-2'))]
      frac_loss = ifelse(length(frac_loss) == 2, (frac_loss[[1]] + frac_loss[[2]]), frac_loss[[1]])
      
      frac_gain = table_out[which(names(table_out) %in% c('1', '2'))]
      frac_gain = ifelse(length(frac_gain) == 2, (frac_gain[[1]] + frac_gain[[2]]), frac_gain[[1]])
      
      out = data.frame(id = id,
                       p.id = p.id,
                       fga = fga,
                       purity = purity,
                       diplogr = diplogr,
                       n_segs = n_segs,
                       frac_below_diplogr = frac_below_diplogr,
                       waterfall_sd = waterfall,
                       qc = qc_facets,
                       frac_loss = (frac_loss) / sum(table(gene$numeric_call)),
                       frac_gain = (frac_gain) / sum(table(gene$numeric_call)),
                       cdkn2a = cdkn2a_call,
                       type = type)
      
    } else {
      file = all_files[grep(pattern = id, x = all_files)]
      file = file[grep(pattern = 'adjusted.*', x = file)]
      filerds = file[grep(pattern = 'rds$', file)]
      data = readRDS(filerds)
      
      ##-- run Tests
      qc = facets_fit_qc(facets_output = data)
      gene = facetsSuite::gene_level_changes(facets_output = data, genome = 'hg19', algorithm = 'em')
      gene$string = paste(qc$wgd, gene$tcn.em, gene$mcn, gene$lcn.em, sep = ':')
      gene = merge(gene, copynumberstates[,c('map_string', 'numeric_call')],  by.x = 'string', by.y = 'map_string', all.x = T)
      cdkn2a_call = gene$numeric_call[which(gene$gene == 'CDKN2A')]
      
      fga = qc$fga
      purity = qc$purity
      diplogr = qc$dipLogR
      n_segs = qc$n_segs
      frac_below_diplogr = qc$frac_below_dipLogR
      waterfall = qc$sd_cnlr_residual
      qc_facets = qc$facets_qc
      
      ##-- alterations
      table_out = table(gene$numeric_call)
      
      frac_loss = table_out[which(names(table_out) %in% c('-1', '-2'))]
      frac_loss = ifelse(length(frac_loss) == 2, (frac_loss[[1]] + frac_loss[[2]]), frac_loss[[1]])
      
      frac_gain = table_out[which(names(table_out) %in% c('1', '2'))]
      frac_gain = ifelse(length(frac_gain) == 2, (frac_gain[[1]] + frac_gain[[2]]), frac_gain[[1]])
      
      out = data.frame(id = id,
                       p.id = p.id,
                       fga = fga,
                       purity = purity,
                       diplogr = diplogr,
                       n_segs = n_segs,
                       frac_below_diplogr = frac_below_diplogr,
                       waterfall_sd = waterfall,
                       qc = qc_facets,
                       frac_loss = (frac_loss) / sum(table(gene$numeric_call)),
                       frac_gain = (frac_gain) / sum(table(gene$numeric_call)),
                       cdkn2a = cdkn2a_call,
                       type = type)
      
    }
    cdkn2a = rbind(cdkn2a, out)
  })
}



##-- CnLR plots
a = readRDS('01_countmatrices/countsMerged____P-0008244-T01-IM5_P-0008244-N01-IM5.dat.gz/P-0008244-T01-IM5_P-0008244-N01-IM5.rds')
b = readRDS('01_countmatrices/s_C_002114_NCAS_dZ_IM5.rg.md.abra.printreads--s_C_002114_S003_d.rg.md.abra.printreads.pileup/s_C_002114_NCAS_dZ_IM5.rg.md.abra.printreads--s_C_002114_S003_d.rds')

a_plot = cnlr_plot(facets_data = a, highlight_gene = 'CDKN2A')
a_plot = a_plot + labs(title = 'P-0008244-T01-IM5')

b_plot = cnlr_plot(facets_data = b, highlight_gene = 'CDKN2A')
b_plot = b_plot + labs(title = 's_C_002114_S003')

dip_cdkn2a = a_plot / b_plot
ggsave_golden(filename = '05_Plots/CDKN2A_diploid.pdf', plot = dip_cdkn2a, width = 10)


##-- look into the coverage distribution:
cov = read.csv('04_MADSEQ/TumorCSF/s_C_002114_S003_d_normed_depth.txt', sep = '\t')
cov = cov[which(cov$seqnames %in% c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7',
                                    'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14',
                                    'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22')), ]

cov$seqnames = factor(cov$seqnames, levels = c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7',
                                               'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14',
                                               'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22'))
csf1 = ggplot(cov, aes(x = normed_depth)) + 
  geom_density(color = "black", fill = "white") +
  scale_x_continuous(breaks = c(0, 100, 200),
                     labels = c(0, 100, 200)) +
  scale_y_continuous(expand = c(0.01, 0)) +
  theme_std() +
  facet_grid(~seqnames) +
  theme(axis.text.x = element_text(angle = 90),
        axis.text.y = element_blank()) +
  labs(y = 'Density') +
  panel_border(size = 0.75, color = 'grey85')
  
tumor = read.csv('04_MADSEQ/TumorCSF/P-0008244-T01-IM5_normed_depth.txt', sep = '\t')
tumor = tumor[which(tumor$seqnames %in% c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7',
                                      'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14',
                                      'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22')), ]
tumor$seqnames = factor(tumor$seqnames, levels = c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7',
                                               'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14',
                                               'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22'))

tumor1 = ggplot(tumor, aes(x = normed_depth)) + 
  geom_density(color = "black", fill = "white") +
  scale_x_continuous(breaks = c(0, 750, 1500),
                     labels = c(0, 750, 1500)) +
  scale_y_continuous(expand = c(0.01, 0)) +
  theme_std() +
  facet_grid(~seqnames) +
  theme(axis.text.x = element_text(angle = 90),
        axis.text.y = element_blank()) +
  labs(y = 'Density') +
  panel_border(size = 0.75, color = 'grey85')


coverage = tumor1 / csf1
ggsave_golden(filename = '05_Plots/CDKN2A_diploid_coverage.pdf', plot = coverage, width = 12)



##----------------+
##-- DeepDeletion
##----------------+
a = readRDS('01_countmatrices/countsMerged____P-0031834-T01-IM6_P-0031834-N01-IM6.dat.gz/P-0031834-T01-IM6_P-0031834-N01-IM6.rds')
b = readRDS('01_countmatrices/s_C_RT2H2R_N901_dZ_IM6.rg.md.abra.printreads__s_C_RT2H2R_S001_d.rg.md.abra.printreads.dat.gz/s_C_RT2H2R_N901_dZ_IM6.rg.md.abra.printreads__s_C_RT2H2R_S001_d.rds')

a_plot = cnlr_plot(facets_data = a, highlight_gene = 'CDKN2A')
a_plot = a_plot + labs(title = 'P-0031834-T01-IM6')

b_plot = cnlr_plot(facets_data = b, highlight_gene = 'CDKN2A')
b_plot = b_plot + labs(title = 's_C_RT2H2R_S001')

deep_cdkn2a = a_plot / b_plot
ggsave_golden(filename = '05_Plots/CDKN2A_deepDeletion.pdf', plot = deep_cdkn2a, width = 10)


##-- Coverage
cov = read.csv('04_MADSEQ/TumorCSF/s_C_RT2H2R_S001_d_normed_depth.txt', sep = '\t')
cov = cov[which(cov$seqnames %in% c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7',
                                    'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14',
                                    'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22')), ]

cov$seqnames = factor(cov$seqnames, levels = c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7',
                                               'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14',
                                               'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22'))
csf1 = ggplot(cov, aes(x = normed_depth)) + 
  geom_density(color = "black", fill = "white") +
  #scale_x_continuous(breaks = c(0, 100, 200),
  #                   labels = c(0, 100, 200)) +
  scale_y_continuous(expand = c(0.01, 0)) +
  theme_std() +
  facet_grid(~seqnames) +
  theme(axis.text.x = element_text(angle = 90),
        axis.text.y = element_blank()) +
  labs(y = 'Density') +
  panel_border(size = 0.75, color = 'grey85')


tumor = read.csv('04_MADSEQ/TumorCSF/P-0031834-T01-IM6_normed_depth.txt', sep = '\t')
tumor = tumor[which(tumor$seqnames %in% c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7',
                                          'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14',
                                          'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22')), ]
tumor$seqnames = factor(tumor$seqnames, levels = c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7',
                                                   'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14',
                                                   'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22'))

tumor1 = ggplot(tumor, aes(x = normed_depth)) + 
  geom_density(color = "black", fill = "white") +
  #scale_x_continuous(breaks = c(0, 750, 1500),
  #                   labels = c(0, 750, 1500)) +
  scale_y_continuous(expand = c(0.01, 0)) +
  theme_std() +
  facet_grid(~seqnames) +
  theme(axis.text.x = element_text(angle = 90),
        axis.text.y = element_blank()) +
  labs(y = 'Density') +
  panel_border(size = 0.75, color = 'grey85')


coverage = tumor1 / csf1
ggsave_golden(filename = '05_Plots/CDKN2A_deepDeletion_coverage.pdf', plot = coverage, width = 12)


##----------------+
## CDKN2A
##----------------+
cdkn2a$type = factor(cdkn2a$type, levels = c('TUMOR', 'CSF'))
waterfall_concordant = ggplot(cdkn2a, aes(x = type, y = waterfall_sd)) +
  geom_boxplot(outlier.size = 0.5) +
  geom_jitter(width = 0.2, size = 1.5) +
  theme_std(base_size = 14, base_line_size = 0.1) +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0, 1.5)) +
  panel_border(size = 2, color = 'black') +
  theme(aspect.ratio = 1.5,
        axis.ticks.y = element_line(linewidth = 0.85)) +
  labs(y = 'Log2 Ratio residual [sd]', x = '', title = 'Failed Sequencing (Waterfall)?')

n_segs_concordant = ggplot(cdkn2a, aes(x = type, y = n_segs)) +
  geom_boxplot(outlier.size = 0.5) +
  geom_jitter(width = 0.2, size = 1.5) +
  theme_std(base_size = 14, base_line_size = 0.1) +
  scale_y_continuous(expand = c(0,0),
                     limits = c(25, 75)) +
  panel_border(size = 2, color = 'black') +
  theme(aspect.ratio = 1.5,
        axis.ticks.y = element_line(linewidth = 0.85)) +
  labs(y = '# Segments', x = '', title = 'Hypersegmentation?')


Purity_concordance = ggplot(cdkn2a, aes(x = type, y = purity)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 1.5) +
  theme_std(base_size = 14, base_line_size = 0.1) +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0, 1),
                     labels = c('0%', '25%', '50%', '75%', '100%')) +
  panel_border(size = 2, color = 'black') +
  theme(aspect.ratio = 1.5,
        axis.ticks.y = element_line(linewidth = 0.85)) +
  labs(y = 'Purity', x = '', title = 'Purity?')


concordant = waterfall_concordant + n_segs_concordant + Purity_concordance
ggsave_golden(filename = '05_Plots/CDKN2A_concordance.pdf', plot = concordant, width = 14)



##----------------+
## missing CDKN2A calls in CSF
##----------------+
discordant_cdkn2a = c('C-001327', 'C-001796', 'C-006881','C-006889',
                      'C-10D913', 'C-2V7CV9', 'C-38YEPA', 'C-5Y6LA4',
                      'C-C4F76K', 'C-H4LNTV')
discordant_cdkn2a = pairing[which(pairing$Patient.ID %in% discordant_cdkn2a & pairing$Paired_Analysis == 'use'), ]

dis_cdkn2a = data.frame()
for(i in 1:nrow(discordant_cdkn2a)){
  try({
    print(i)
    id = discordant_cdkn2a$Sample.ID[i]
    fit = discordant_cdkn2a$Fit[i]
    p.id = unique(discordant_cdkn2a$Patient.ID[i])
    type = unique(discordant_cdkn2a$TYPE[i])
    
    if(fit == 'full'){
      file = all_files[grep(pattern = id, x = all_files)]
      filerds = file[grep(pattern = '.rds$', x = file)]
      data = readRDS(filerds)
      
      ##-- run Tests
      qc = facets_fit_qc(facets_output = data)
      gene = facetsSuite::gene_level_changes(facets_output = data, genome = 'hg19', algorithm = 'em')
      gene$string = paste(qc$wgd, gene$tcn.em, gene$mcn, gene$lcn.em, sep = ':')
      gene = merge(gene, copynumberstates[,c('map_string', 'numeric_call')],  by.x = 'string', by.y = 'map_string', all.x = T)
      cdkn2a_call = gene$numeric_call[which(gene$gene == 'CDKN2A')]
      
      fga = qc$fga
      purity = qc$purity
      diplogr = qc$dipLogR
      n_segs = qc$n_segs
      frac_below_diplogr = qc$frac_below_dipLogR
      waterfall = qc$sd_cnlr_residual
      qc_facets = qc$facets_qc
      
      ##-- alterations
      table_out = table(gene$numeric_call)
      
      frac_loss = table_out[which(names(table_out) %in% c('-1', '-2'))]
      frac_loss = ifelse(length(frac_loss) == 2, (frac_loss[[1]] + frac_loss[[2]]), frac_loss[[1]])
      
      frac_gain = table_out[which(names(table_out) %in% c('1', '2'))]
      frac_gain = ifelse(length(frac_gain) == 2, (frac_gain[[1]] + frac_gain[[2]]), frac_gain[[1]])
      
      out = data.frame(id = id,
                       p.id = p.id,
                       fga = fga,
                       purity = purity,
                       diplogr = diplogr,
                       n_segs = n_segs,
                       frac_below_diplogr = frac_below_diplogr,
                       waterfall_sd = waterfall,
                       qc = qc_facets,
                       frac_loss = (frac_loss) / sum(table(gene$numeric_call)),
                       frac_gain = (frac_gain) / sum(table(gene$numeric_call)),
                       cdkn2a = cdkn2a_call,
                       type = type)
      
    } else {
      file = all_files[grep(pattern = id, x = all_files)]
      file = file[grep(pattern = 'adjusted.*', x = file)]
      filerds = file[grep(pattern = 'rds$', file)]
      data = readRDS(filerds)
      
      ##-- run Tests
      qc = facets_fit_qc(facets_output = data)
      gene = facetsSuite::gene_level_changes(facets_output = data, genome = 'hg19', algorithm = 'em')
      gene$string = paste(qc$wgd, gene$tcn.em, gene$mcn, gene$lcn.em, sep = ':')
      gene = merge(gene, copynumberstates[,c('map_string', 'numeric_call')],  by.x = 'string', by.y = 'map_string', all.x = T)
      cdkn2a_call = gene$numeric_call[which(gene$gene == 'CDKN2A')]
      
      fga = qc$fga
      purity = qc$purity
      diplogr = qc$dipLogR
      n_segs = qc$n_segs
      frac_below_diplogr = qc$frac_below_dipLogR
      waterfall = qc$sd_cnlr_residual
      qc_facets = qc$facets_qc
      
      ##-- alterations
      table_out = table(gene$numeric_call)
      
      frac_loss = table_out[which(names(table_out) %in% c('-1', '-2'))]
      frac_loss = ifelse(length(frac_loss) == 2, (frac_loss[[1]] + frac_loss[[2]]), frac_loss[[1]])
      
      frac_gain = table_out[which(names(table_out) %in% c('1', '2'))]
      frac_gain = ifelse(length(frac_gain) == 2, (frac_gain[[1]] + frac_gain[[2]]), frac_gain[[1]])
      
      out = data.frame(id = id,
                       p.id = p.id,
                       fga = fga,
                       purity = purity,
                       diplogr = diplogr,
                       n_segs = n_segs,
                       frac_below_diplogr = frac_below_diplogr,
                       waterfall_sd = waterfall,
                       qc = qc_facets,
                       frac_loss = (frac_loss) / sum(table(gene$numeric_call)),
                       frac_gain = (frac_gain) / sum(table(gene$numeric_call)),
                       cdkn2a = cdkn2a_call,
                       type = type)
      
    }
    dis_cdkn2a = rbind(dis_cdkn2a, out)
  })
}

cdkn2a_start = 21967751
cdkn2a_end = 21995323

bb = aa[which(aa$chrom == 9 & aa$maploc >= cdkn2a_start & aa$maploc <= cdkn2a_end), ]
plot(density(bb$cnlr))

boxplot(dis_cdkn2a$waterfall_sd ~ dis_cdkn2a$type)

cnlr_plot(facets_data = data)

aa = data$snps
head(aa)







egfr
C-001327
C-006889
C-006891
C-38YEPA
C-5Y6LA4




dfa
## Cross-check geneLevel CSF calls with cBIO and prepare Plot
## 
## start: 09/13/2022
## chris-Kreitzer
## 

clean()
gc()
setup(working.path = '~/Documents/MSKCC/Subhi/CSF/')
source('~/Documents/GitHub/CSF/Scripts/UtilityFunctions.R')
source('~/Documents/GitHub/CSF/Scripts/plot_theme.R')


cBIO = read.csv(file = '~/Documents/MSKCC/Subhi/CSF/Data/cBIO_template.txt', sep = '\t')
gene_level_path = 'C_DA3H4U/gene_level_out.txt'
gene_level_data = read.csv(gene_level_path, sep = '\t')


#' add cBIO annotation; 
#' manual change of genes
cBIO$group = 'cBIO-DMP1'
cBIO$copy_state[which(cBIO$gene == 'CDK4')] = 'GAIN'
cBIO$copy_state[which(cBIO$gene == 'MDM2')] = 'GAIN'
cBIO$copy_state[!cBIO$gene %in% c('CDK4', 'MDM2')] = 'DIPLOID'

cBIO2 = cBIO
cBIO2$group = 'cBIO-DMP2'
cBIO2$copy_state = 'DIPLOID'


gene_level_data = rbind(cBIO, cBIO2, gene_level_data)
gene_level_data$group[which(gene_level_data$group == 1)] = 'Facets-DMP1'
gene_level_data$group[which(gene_level_data$group == 2)] = 'Facets-DMP2'
gene_level_data$group[which(gene_level_data$group == 3)] = 'CSF1'
gene_level_data$group[which(gene_level_data$group == 4)] = 'CSF2'
gene_level_data$group = factor(gene_level_data$group, levels = c('cBIO-DMP1', 'Facets-DMP1', 'cBIO-DMP2', 'Facets-DMP2', 'CSF1', 'CSF2'))
gene_level_data$clonality = ifelse(gene_level_data$clonality == 'clonal', '', '*')
gene_level_data$copy_state[which(gene_level_data$copy_state == 'CNLOH')] = 'LOSS'
gene_level_data$copy_state[which(gene_level_data$copy_state == 'HOMDEL')] = 'LOSS'
gene_level_data$copy_state[which(gene_level_data$copy_state == 'AMP')] = 'GAIN'
gene_level_data$copy_state[which(gene_level_data$copy_state == 'HETLOSS')] = 'LOSS'
gene_level_data$copy_state[which(gene_level_data$copy_state == 'GAIN (many states)')] = 'GAIN'
gene_level_data$copy_state[which(gene_level_data$copy_state == 'CNLOH & GAIN')] = 'GAIN'
gene_level_data$copy_state[which(gene_level_data$copy_state == 'AMP (many states)')] = 'GAIN'
gene_level_data$copy_state[which(gene_level_data$copy_state == 'TETRAPLOID')] = 'GAIN'


DMP_C_DA3 = ggplot(gene_level_data, aes(x = group, y = gene, fill = copy_state)) +
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

DMP_C_DA3
ggsave(filename = 'C_DA3H4U/geneLevel_DMPnormal.pdf', plot = DMP_C_DA3, device = 'pdf', width = 4, height = 8)

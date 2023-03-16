##----------------+
## CnLR calls for disease
## defining genes in CSF and Tumor
##----------------+
##
## start: 03/16/2023
## chris-kreitzer


clean()
setwd('~/Documents/MSKCC/11_CSF/')
source('~/Documents/GitHub/Y_chromosome_loss/PanCancer/Scripts/UtilityFunctions.R')
all_files = list.files('01_countmatrices/', pattern = 'rds$', recursive = T, full.names = T)
GOIs = c('CDKN2A','CDK4','CDK6','PTEN','EGFR','RET',
         'GATA3','PDGFRA','KIT','KDR','MET','MDM2',
         'MDM4','RB1','TP53')

database = readxl::read_excel('00_Data/Database_Final_march12.xlsx')
database = database[!database$Fit %in% c('N/Avail'), ]


##-- RUNNING CSF cohort
csf_cohort = as.data.frame(database[which(database$TYPE == 'CSF'), c('Sample.ID', 'Fit', 'TYPE')])
csf_out = data.frame()
for(i in 1:nrow(csf_cohort)){
  try({
    print(i)
    file = all_files[grep(pattern = csf_cohort$Sample.ID[i], x = all_files)]
    file = ifelse(length(file) == 1, file[1], 
                  ifelse(length(file) > 1 & any(grepl(pattern = 'adjusted', x = file)), 
                         file[grep(pattern = 'adjusted', x = file)], file[1]))
    datain = readRDS(file = file)
    
    diplogr = datain$dipLogR
    segs = datain$segs
    
    gene = facetsSuite::gene_level_changes(facets_output = datain, genome = 'hg19', algorithm = 'em')
    gene = gene[which(gene$gene %in% GOIs), c('gene', 'chrom', 'seg', "median_cnlr_seg", 'cn_state')]
    
    if(nrow(gene) == 0) next
    
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
    gene_new$id = csf_cohort$Sample.ID[i]
    gene_new$Fit = csf_cohort$Fit[i]
    gene_new$TYPE = csf_cohort$TYPE[i]
    
    csf_out = rbind(csf_out, gene_new)
    rm(segs, gene, datain, file, gene_new, diplogr, cnlr_distribution)
  })
}
csf_out = csf_out[,c(6,8,7,1,3,5)]


##-- Tumor Cohort
tumor_cohort = as.data.frame(database[which(database$TYPE == 'TUMOR'), c('Sample.ID', 'Fit', 'TYPE')])
tumor_out = data.frame()
for(i in 1:nrow(tumor_cohort)){
  try({
    print(i)
    file = all_files[grep(pattern = tumor_cohort$Sample.ID[i], x = all_files)]
    file = ifelse(length(file) == 1, file[1], 
                  ifelse(length(file) > 1 & any(grepl(pattern = 'adjusted', x = file)), 
                         file[grep(pattern = 'adjusted', x = file)], file[1]))
    datain = readRDS(file = file)
    
    diplogr = datain$dipLogR
    segs = datain$segs
    
    gene = facetsSuite::gene_level_changes(facets_output = datain, genome = 'hg19', algorithm = 'em')
    gene = gene[which(gene$gene %in% GOIs), c('gene', 'chrom', 'seg', "median_cnlr_seg", 'cn_state')]
    
    if(nrow(gene) == 0) next
    
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
    gene_new$id = tumor_cohort$Sample.ID[i]
    gene_new$Fit = tumor_cohort$Fit[i]
    gene_new$TYPE = tumor_cohort$TYPE[i]
    
    tumor_out = rbind(tumor_out, gene_new)
    rm(segs, gene, datain, file, gene_new, diplogr, cnlr_distribution)
  })
}

tumor_out = tumor_out[,c(6,8,7,1,3,5)]



CnLR_calls = rbind(csf_out, tumor_out)
write.table(x = CnLR_calls, '00_Data/CnLR_calls_all_03162023.txt', sep = '\t', row.names = F, quote = F)




##----------------+
## Make binary matrix of calls
##----------------+
##----------------+
## transform numeric 
## alteration matrix in 
## binary format; no FILTERS considered
##----------------+
alterations = read.csv('00_Data/CnLR_calls_all_03162023.txt', sep = '\t')

## modify matrix
alteration_matrix = setNames(data.frame(matrix(ncol = 15, nrow = 0)), unique(alterations$gene))

all_out = data.frame()
counter = 1
for(id in unique(alterations$id)){
  data.mut.sub = alterations[which(alterations$id == id), ]
  if(nrow(data.mut.sub) != 0){
    for(j in unique(data.mut.sub$gene)){
      if(j %in% colnames(alteration_matrix)){
        alteration_matrix[counter, j] = data.mut.sub$call[which(data.mut.sub$gene == j)][1]
      }
    }
    alteration_matrix[counter, 'Sample.ID'] = id
    counter = counter + 1
    all_out = rbind(all_out, alteration_matrix)
  }
}

all_out = all_out[!duplicated(all_out), ]
all_out = t(all_out)
colnames(all_out) = all_out[nrow(all_out), ]
all_out = all_out[-nrow(all_out), ]

write.table(x = all_out, file = '~/Documents/MSKCC/11_CSF/00_Data/Tumor_CSF_binary_noFilter.txt', sep = '\t')


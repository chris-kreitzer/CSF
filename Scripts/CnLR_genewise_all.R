##----------------+
## Gene assignments:
## binary, cnlr_segment
## and cnlr_gene
##----------------+

clean()
gc()
.rs.restartR()
setwd(dir = '~/Documents/MSKCC/Subhi/CSF/')

GOI = genes_hg19[which(genes_hg19$gene %in% c('CDKN2A',
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
                                              'TP53', 
                                              'ATRX', 
                                              'BRAF',
                                              'PIK3CA', 
                                              'PIK3R1',
                                              'PIK3R2')), ]

folders = list.files(path = '.')
folders = folders[grepl(pattern = 'C-', x = folders)]

all_out = data.frame()
for(i in unique(folders)){
  try({
    sub_dirs = list.dirs(path = i, full.names = T, recursive = F)
    if(length(sub_dirs) == 0) next
    else {
      for(j in unique(sub_dirs)){
        Rdata = list.files(pattern = '.Rdata$', path = paste0(j, '/'), full.names = T)
        load(file = Rdata)
        
        #' compile the whole FACETS output
        name = basename(j)
        print(name)
        
        facets_out = list(
          snps = out$jointseg,
          segs = fit$cncf,
          purity = as.numeric(fit$purity),
          ploidy = as.numeric(fit$ploidy),
          dipLogR = out$dipLogR,
          alBalLogR = out$alBalLogR,
          flags = out$flags,
          em_flags = fit$emflags,
          loglik = fit$loglik)
        
        gene_level = facetsSuite::gene_level_changes(facets_output = facets_out, genome = 'hg19')
        
        for(z in 1:nrow(GOI)){
          gene = GOI$gene[z]
          cnlr = out$jointseg$cnlr[which(out$jointseg$chrom == GOI$chrom[z] & 
                                           out$jointseg$maploc >= GOI$start[z] & 
                                           out$jointseg$maploc <= GOI$end[z])]
          cnlr_mean = mean(cnlr, na.rm = T)
          cnlr_segment = gene_level$median_cnlr_seg[which(gene_level$gene == GOI$gene[z])]
          gene_cnlr = data.frame(id = name,
                                 gene = gene,
                                 cnlr_gene = cnlr_mean,
                                 cnlr_segment = cnlr_segment)
          
          all_out = rbind(all_out, gene_cnlr)
        }
      }
    }
  })
}

# write.table(all_out, file = 'Data/Final/CnLR_genewise_all.txt', sep = '\t')

binary_alterations = read.csv('Data/Final/CSF_binary_Filtered_QC_True.txt', sep = '\t')
colnames(binary_alterations) = gsub(pattern = '\\.', replacement = '-', colnames(binary_alterations))
samples_pass = unique(colnames(binary_alterations))
all_out = all_out[which(all_out$id %in% samples_pass), ]
all_out$Facets_call = NA

for(i in 1:nrow(all_out)){
  gg = binary_alterations[which(colnames(binary_alterations) == all_out$id[i])]
  ggg = gg[which(row.names(gg) == all_out$gene[i]), ]
  all_out$Facets_call[i] = ggg
}


all_out$Facets_call[which(all_out$Facets_call == -1)] = 'Loss'
all_out$Facets_call[which(all_out$Facets_call == -2)] = 'Homozygous Deletion'
all_out$Facets_call[which(all_out$Facets_call == 1)] = 'Gain'
all_out$Facets_call[which(all_out$Facets_call == 2)] = 'Amplification'
all_out$Facets_call[which(all_out$Facets_call == 0)] = 'Diploid'

write.table(all_out, file = 'Data/Final/CnLR_genewise_all.txt', sep = '\t')





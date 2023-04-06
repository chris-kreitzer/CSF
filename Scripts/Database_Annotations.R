##----------------+
## Make annotations for the final sheed
## - Nic Socci's calls
## - imitated MSK-pipeline
## - Facets
##----------------+
## 
## start: 04/06/2023
## chris-kreitzer


setwd('~/Documents/MSKCC/11_CSF/')
GOIs = c('CDKN2A','CDK4','CDK6','PTEN','EGFR','PDGFRA','KIT','KDR','MET','MDM2','MDM4','RB1','NF1','TP53','FGF4', 'FGF19')
database = read.csv('00_Data/Chris_Subhi_Tab_fixed_USE.txt', sep = '\t')
database_old = readxl::read_excel('00_Data/Chris_Subhi_Data_April_5_2023.xlsx')
database_old = database_old[,c('Sample.ID', 'FacetsCountfile', 'Fit', 'TumorBam', 'NormalBam', "arm_7p", 'arm_7q', "arm_10p", "arm_10q", 'purity', "ploidy", 'wgd', 'fga', "SCNA_Chris")]


##-- Nic calls
a = readxl::read_excel('06_Nic_Socci_r_001/seqCNA/IM-gbm_IM5/COHORT/IM-gbm_IM5___GeneTable_Vx_1_1.xlsx')
b = readxl::read_excel('06_Nic_Socci_r_001/seqCNA/IM-gbm_IM6/COHORT/IM-gbm_IM6___GeneTable_Vx_1_1.xlsx')
c = readxl::read_excel('06_Nic_Socci_r_001/seqCNA/IM-gbm_IM7/COHORT/IM-gbm_IM7___GeneTable_Vx_1_1.xlsx')

Nic = rbind(a,b,c)
Nic$call = ifelse(Nic$seg.mean <= -2, -2,
                  ifelse(Nic$seg.mean >= 2, 2, 0))

Nic = Nic[which(Nic$gene %in% GOIs), ]
write.table(x = Nic, file = '06_Nic_Socci_r_001/Nic_comprehensive_CNA_calls.txt', sep = '\t', row.names = F, quote = F)

##-- add to data table:
Nic_com = data.frame()
for(i in unique(Nic$Tumor)){
  genes = paste(Nic$gene[which(Nic$Tumor == i & Nic$call != 0)], collapse = ',')
  out = data.frame(id = i,
                   SCNA_Nic = genes)
  Nic_com = rbind(Nic_com, out)
}

Nic_com 


##-- Imitated MSK-pipeline
Impact = read.csv('00_Data/MSK_like_CNA_calls.txt', sep = '\t')
Impact = Impact[which(Impact$flag_pileup == FALSE), ]

chris_com = data.frame()
for(i in unique(Impact$id)){
  genes = paste(Impact$gene[which(Impact$id == i & Impact$call != 0)], collapse = ',')
  out = data.frame(id = i,
                   SCNA_chris = genes)
  chris_com = rbind(chris_com, out)
}



##-- Merge with original Table (04/06/2023)
data = merge(database, database_old, by.x = 'Sample.ID', by.y = 'Sample.ID', all.x = T)
data = merge(data, Nic_com, by.x = 'Sample.ID', by.y = 'id', all.x = T)
data = merge(data, chris_com, by.x = 'Sample.ID', by.y = 'id', all.x = T)

write.table(x = data, file = '00_Data/Chris_Subhi_Tab_fixed_USE_04062023.txt', sep = '\t', row.names = F, quote = F)


#-----------------+
# Make binary table 
# for different datasets
# - starting with Nic's calls
#-----------------+
alterations = read.csv('06_Nic_Socci_r_001/Nic_comprehensive_CNA_calls.txt', sep = '\t')

## modify matrix
alteration_matrix = setNames(data.frame(matrix(ncol = length(unique(alterations$gene)), nrow = 0)), unique(alterations$gene))

all_out = data.frame()
counter = 1
for(id in unique(alterations$Tumor)){
  data.mut.sub = alterations[which(alterations$Tumor == id), ]
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
all_out[is.na(all_out)] = 0
write.table(x = all_out, file = '~/Documents/MSKCC/11_CSF/00_Data/Nic_Socci_binary.txt', sep = '\t')



##-- chris imitated pipeline
alterations = read.csv('00_Data/MSK_like_CNA_calls.txt', sep = '\t')
alterations = alterations[which(alterations$flag_pileup == FALSE), ]

## modify matrix
alteration_matrix = setNames(data.frame(matrix(ncol = length(unique(alterations$gene)), nrow = 0)), unique(alterations$gene))

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
all_out[is.na(all_out)] = 0
write.table(x = all_out, file = '~/Documents/MSKCC/11_CSF/00_Data/chris_imitated_MSK_binary.txt', sep = '\t')






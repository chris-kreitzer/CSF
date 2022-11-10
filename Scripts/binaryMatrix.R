##---------------------------
## create binary Matrix for 
## CSF alterations:
##---------------------------
##
## start: 10/07/2022
## update: 10/12/2022
## update: 10/27/2022
## update: 11/10/2022
## chris-kreitzer


clean()
gc()
.rs.restartR()
setup(working.path = '~/Documents/MSKCC/Subhi/CSF/')
folders = list.files(path = '.')
folders = folders[grepl(pattern = 'C-', x = folders)]

##' Definitions
copynumberstates = facetsSuite:::copy_number_states
GOI = c("CDKN2A", "CDKN2B", "MTAP", 'EGFR', 'CDK4', 'PDGFRA', 'PTEN', 
        'KIT', 'MDM2', 'KDR', 'MDM4', 'RB1', 'MET', 'NF1', 'CDK6', 
        'TP53', 'KRAS', 'ATRX', 'FGF3', 'FGF4', 'FGF19')
gene_filter_states = c('suppress_segment_too_large', 
                       'suppress_likely_unfocal_large_gain', 
                       'suppress_large_homdel')


alterations = data.frame()
IGV_all = data.frame()
for(i in unique(folders)){
  try({
    sub_dirs = list.dirs(path = i, full.names = T, recursive = F)
    if(length(sub_dirs) == 0) next
    else {
      for(j in unique(sub_dirs)){
        Rdata = list.files(pattern = '.Rdata$', path = paste0(j, '/'), full.names = T)
        load(file = Rdata)
        fit$cncf$lcn[fit$cncf$tcn == 1] = 0
        
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
        
        #' gene_level_changes() and other metrics 
        gene_out = facetsSuite::gene_level_changes(facets_output = facets_out, genome = 'hg19', algorithm = 'cncf')
        IGV = facetsSuite::format_igv_seg(facets_output = facets_out, sample_id = j, normalize = T)
        gene_out = gene_out[which(gene_out$gene %in% GOI), ]
        fcna_output = facetsSuite::calculate_fraction_cna(facets_out$segs, facets_out$ploidy, genome = 'hg19', algorithm = 'cncf')
        wgd = fcna_output$genome_doubled
        
        #---------+
        # loop through genes
        #---------+
        for(k in 1:nrow(gene_out)){
          filter = gene_out$filter[k]
          gene = gene_out$gene[k]
          call = gene_out$cn_state[k]
          numeric_call = ifelse(call == 'GAIN (many states)' & wgd, -1,
                                ifelse(call == 'TETRAPLOID' & wgd, 0, 
                                       unique(copy_number_states$numeric_call[which(copy_number_states$call == call)])))
          
          gene_final = data.frame(id = name,
                                  gene = gene,
                                  call = call,
                                  n_call = numeric_call)
          
          alterations = rbind(alterations, gene_final)
        }
      }
      IGV_all = rbind(IGV_all, IGV)
      rm(fit, out, name, facets_out, gene_out, gene_final)
    }
  })
}









##-----------------
## work on ATRX:
##-----------------
folders = list.files(path = '.')
folders = folders[grepl(pattern = 'C-', x = folders)]
sample_clinical = readxl::read_excel('Data/FINAL_samples/CSF_Lastest_07102022.xlsx')
sample_clinical_short = sample_clinical[, c('Patient ID', 'Sample ID', 'Sex')]
sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-09T7CL')] = 'M'
sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-001615')] = 'M'
sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-4173R4')] = 'F'
sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-56D79J')] = 'M'
sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-5XD5TE')] = 'M'
sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-65DNPK')] = 'M'
sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-6JMECR')] = 'M'
sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-9E8T53')] = 'M'
sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-A32F9N')] = 'M'
sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-NRMAAD')] = 'M'
sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-UMDEWL')] = 'M'
sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-DVCFA2')] = 'F'
sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-D2P3RM')] = 'M'
sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-PTK9RU')] = 'M'
sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-NC3N8D')] = 'F'
sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-D7MUWF')] = 'M'
colnames(sample_clinical_short) = c('PatientID', 'SampleID', 'Sex')

all_out = data.frame()
for(i in unique(folders)){
  try({
    sub_dirs = list.dirs(path = i, full.names = T, recursive = F)
    if(length(sub_dirs) == 0) next
    else {
      for(j in unique(sub_dirs)){
        Rdata = list.files(pattern = '.Rdata$', path = paste0(j, '/'), full.names = T)
        load(file = Rdata)
        fit$cncf$cf = NULL
        fit$cncf$tcn = NULL
        fit$cncf$lcn = NULL
        fit$cncf = cbind(fit$cncf, cf = out$out$cf, tcn = out$out$tcn, lcn = out$out$lcn)
        fit$cncf$lcn[fit$cncf$tcn == 1] = 0
        fit$cncf$lcn.em[fit$cncf$tcn.em == 1] = 0
        
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
        
        gene_out = facetsSuite::gene_level_changes(facets_output = facets_out,
                                                   genome = 'hg19')
        ATRX = gene_out[which(gene_out$gene == 'ATRX'), ]
        wgd = facets_fit_qc(facets_output = facets_out)$wgd
        
        
        ##-----------
        if(wgd & sample_clinical_short$Sex[which(sample_clinical_short$SampleID == name)] == 'F'){
          copyStates = copynumberstates[which(copynumberstates$wgd == T), ]
          tcn = ATRX$tcn.em
          lcn = ATRX$lcn.em
          filter = ATRX$filter
          call = ifelse(filter %in% gene_filter_states, NA,
                        ifelse(!filter %in% gene_filter_states & is.na(lcn), 
                               copyStates$numeric_call[which(copyStates$tcn == tcn &
                                                               is.na(copyStates$lcn))],
                               ifelse(!filter %in% gene_filter_states & !is.na(lcn),
                                      copyStates$numeric_call[which(copyStates$tcn == tcn &
                                                                      copyStates$lcn == lcn)], NA)))
          out = data.frame(id = name,
                           gene = 'ATRX', 
                           call = call)
          
        } else if (!wgd & sample_clinical_short$Sex[which(sample_clinical_short$SampleID == name)] == 'F'){
          copyStates = copynumberstates[which(copynumberstates$wgd == F), ]
          tcn = ATRX$tcn.em
          lcn = ATRX$lcn.em
          filter = ATRX$filter
          call = ifelse(filter %in% gene_filter_states, NA,
                        ifelse(!filter %in% gene_filter_states & is.na(lcn), 
                               copyStates$numeric_call[which(copyStates$tcn == tcn &
                                                               is.na(copyStates$lcn))],
                               ifelse(!filter %in% gene_filter_states & !is.na(lcn),
                                      copyStates$numeric_call[which(copyStates$tcn == tcn &
                                                                      copyStates$lcn == lcn)], NA)))
          out = data.frame(id = name,
                           gene = 'ATRX', 
                           call = call)
          
        } else if (wgd & sample_clinical_short$Sex[which(sample_clinical_short$SampleID == name)] == 'M'){
          tcn = ATRX$tcn.em
          lcn = ATRX$lcn.em
          filter = ATRX$filter
          call = ifelse(filter %in% gene_filter_states, NA,
                        ifelse(!filter %in% gene_filter_states & lcn == 0 & tcn == 1, -1,
                               ifelse(!filter %in% gene_filter_states & is.na(lcn) & tcn == 1, -1,
                                      ifelse(!filter %in% gene_filter_states & lcn == 0 & tcn > 1, 0,
                                             ifelse(!filter %in% gene_filter_states & is.na(lcn) & tcn > 1, 0,
                                                    ifelse(!filter %in% gene_filter_states & lcn == 0 & tcn == 0, -1, 
                                                           ifelse(!filter %in% gene_filter_states & is.na(lcn) & tcn > 2, 1, NA)))))))
                                             
                               
                               
          out = data.frame(id = name,
                           gene = 'ATRX', 
                           call = call)
          
        } else if (!wgd & sample_clinical_short$Sex[which(sample_clinical_short$SampleID == name)] == 'M'){
          tcn = ATRX$tcn.em
          lcn = ATRX$lcn.em
          filter = ATRX$filter
          call = ifelse(filter %in% gene_filter_states, NA,
                        ifelse(!filter %in% gene_filter_states & is.na(lcn) & tcn > 1, 1,
                               ifelse(!filter %in% gene_filter_states & is.na(lcn) & tcn == 1, 0,
                                      ifelse(!filter %in% gene_filter_states & lcn == 0 & tcn > 1, 1,
                                             ifelse(!filter %in% gene_filter_states & lcn == 0 & tcn == 1, 0,
                                                    ifelse(!filter %in% gene_filter_states & is.na(lcn) & tcn == 0, -1,
                                                           ifelse(!filter %in% gene_filter_states & lcn == 0 & tcn == 0, -1, NA)))))))
                                             
          
          out = data.frame(id = name,
                           gene = 'ATRX', 
                           call = call)
                               
        }
        
        all_out = rbind(all_out, out)
      }
    }
  })
}


##-----------------
## MERGE alterations (general)
## and ATRX
##-----------------
alterations = alterations[which(alterations$gene != 'ATRX'), ]
alterations_all = rbind(alterations, all_out)

## modify matrix
alteration_matrix = setNames(data.frame(matrix(ncol = 21,
                                               nrow = 0)), 
                                 unique(alterations_all$gene))

all_out = data.frame()
counter = 1

for(id in unique(alterations_all$id)){
  data.mut.sub = alterations_all[which(alterations_all$id == id), ]
  if(nrow(data.mut.sub) != 0){
    for(j in unique(data.mut.sub$gene)){
      if(j %in% colnames(alteration_matrix)){
        alteration_matrix[counter, j] = data.mut.sub$call[which(data.mut.sub$gene == j)][1]
      }
    }
    
  }
  alteration_matrix[counter, 'Sample.ID'] = id
  counter = counter + 1
  all_out = rbind(all_out, alteration_matrix)
}
  
all_out = all_out[!duplicated(all_out), ]
all_out = t(all_out)
colnames(all_out) = all_out[nrow(all_out), ]
all_out = all_out[-nrow(all_out), ]

write.table(x = all_out, file = '~/Documents/MSKCC/Subhi/CSF/Data/FINAL_samples/CSF_binary_all.txt', sep = '\t')



##-----------------
## Format IGV:
##-----------------
IGV_all = data.frame()
for(i in unique(folders)){
  try({
    sub_dirs = list.dirs(path = i, full.names = T, recursive = F)
    if(length(sub_dirs) == 0) next
    else {
      for(j in unique(sub_dirs)){
        Rdata = list.files(pattern = '.Rdata$', path = paste0(j, '/'), full.names = T)
        load(file = Rdata)
        fit$cncf$cf = NULL
        fit$cncf$tcn = NULL
        fit$cncf$lcn = NULL
        fit$cncf = cbind(fit$cncf, cf = out$out$cf, tcn = out$out$tcn, lcn = out$out$lcn)
        fit$cncf$lcn[fit$cncf$tcn == 1] = 0
        fit$cncf$lcn.em[fit$cncf$tcn.em == 1] = 0
        
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
        
        IGV = facetsSuite::format_igv_seg(facets_output = facets_out, sample_id = j, normalize = T)
        IGV_all = rbind(IGV_all, IGV)
      }
    }
  })
}

IGV_all$ID = basename(IGV_all$ID)
write.table(x = IGV_all, file = '~/Documents/MSKCC/Subhi/CSF/Data/FINAL_samples/IGV_all.seg', sep = '\t', row.names = F)        



##-----------------
## Filter for QC true
## samples
##-----------------
sample_match = read.csv('Data/FINAL_samples/sample_match.txt', sep = '\t')
samples_pass = sample_match[which(sample_match$fit == 'pass'), ]
samples_pass = samples_pass[!grepl(pattern = 'N01', samples_pass$sample), ]
alterations_pass = all_out[,which(colnames(all_out) %in% samples_pass$sample)]
IGV_pass = IGV_all[which(IGV_all$ID %in% samples_pass$sample), ]

write.table(alterations_pass, file = 'Data/FINAL_samples/CSF_binary_QC_true.txt', sep = '\t')
write.table(IGV_pass, file = 'Data/FINAL_samples/IGV_QC_true.seg', sep = '\t', row.names = F)










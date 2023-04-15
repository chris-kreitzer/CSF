##---------------------------
## create binary Matrix for 
## CSF alterations:
##---------------------------
##
## start: 10/07/2022
## update: 10/12/2022
## update: 10/27/2022
## update: 11/10/2022
## update: 11/13/2022
## chris-kreitzer


clean()
gc()
.rs.restartR()
setup(working.path = '~/Documents/MSKCC/Subhi/CSF/')
folders = list.files(path = '.')
folders = folders[grepl(pattern = 'C-', x = folders)]
Sex_annotation = read.csv('Data/Final/Samples_Sex_Annotation.txt', sep = '\t')
Sex_annotation$path = paste(Sex_annotation$PatientID, Sex_annotation$SampleID, sep = '/')

##' Definitions
copynumberstates = facetsSuite:::copy_number_states
GOI = c("CDKN2A", 'EGFR', 'CDK4', 'PDGFRA', 'PTEN', 
        'KIT', 'MDM2', 'KDR', 'MDM4', 'RB1', 'MET', 
        'NF1', 'CDK6', 'TP53', 'ATRX', 'PIK3CA', 
        'BRAF', 'PIK3R1', 'PIK3R2')

gene_filter_states = c('suppress_segment_too_large', 
                       'suppress_likely_unfocal_large_gain', 
                       'suppress_large_homdel')


##----------------+
## loop through all 
## folders; caution! SEX
##----------------+
alterations = data.frame()
IGV_all = data.frame()
for(i in unique(folders)){
  try({
    sub_dirs = list.dirs(path = i, full.names = T, recursive = F)
    if(length(sub_dirs) == 0) next
    else {
      sub_dirs = sub_dirs[which(sub_dirs %in% Sex_annotation$path)]
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
        gene_out$cn_state = ifelse(gene_out$cn_state == 'INDETERMINATE' & gene_out$tcn > 4, 'AMP', gene_out$cn_state)
        fcna_output = facetsSuite::calculate_fraction_cna(facets_out$segs, facets_out$ploidy, genome = 'hg19', algorithm = 'cncf')
        wgd = fcna_output$genome_doubled
        
        #---------+
        # loop through genes
        #---------+
        for(k in 1:nrow(gene_out)){
          filter = gene_out$filter[k]
          gene = gene_out$gene[k]
          if(gene == 'ATRX'){
            if(Sex_annotation$Sex[which(Sex_annotation$SampleID == name)] == 'M'){
              call = gene_out$cn_state[k]
              numeric_call = ifelse(!wgd & call %in% c('HETLOSS'), 0,
                                    ifelse(!wgd & call %in% c('HOMDEL'), -2,
                                           ifelse(!wgd & call %in% c('CNLOH', 'DIPLOID or CNLOH', 'CNLOH & GAIN', 'GAIN (many states)',
                                                                     'CNLOH & GAIN'), 1, 
                                                  ifelse(wgd & call %in% c('LOSS BEFORE', 'LOSS BEFORE or DOUBLE LOSS AFTER'), 0,
                                                         ifelse(wgd & call %in% c('LOSS BEFORE & AFTER'), -1,
                                                                ifelse(wgd & call %in% c('HOMDEL'), -2,
                                                                       ifelse(wgd & call %in% c('CNLOH BEFORE & LOSS', 'LOSS (many states)',
                                                                                                'CNLOH BEFORE', 'TETRAPLOID or CNLOH BEFORE'), 1, 2)))))))
            } else {
              call = gene_out$cn_state[k]
              numeric_call = ifelse(call == 'GAIN (many states)' & wgd, -1,
                                    ifelse(call == 'TETRAPLOID' & wgd, 0, 
                                           unique(copynumberstates$numeric_call[which(copynumberstates$call == call)])))
            }
            
            gene_final = data.frame(id = name,
                                    gene = gene,
                                    filter = filter,
                                    call = call,
                                    n_call = numeric_call)
            
          } else {
            call = gene_out$cn_state[k]
            numeric_call = ifelse(call == 'GAIN (many states)' & wgd, -1,
                                  ifelse(call == 'TETRAPLOID' & wgd, 0, 
                                         unique(copynumberstates$numeric_call[which(copynumberstates$call == call)])))
            
            gene_final = data.frame(id = name,
                                    gene = gene,
                                    filter = filter,
                                    call = call,
                                    n_call = numeric_call)
          }
          alterations = rbind(alterations, gene_final)
        }
        IGV_all = rbind(IGV_all, IGV)
        rm(fit, out, name, facets_out, gene_out, gene_final)
      }
    }
  })
}

IGV_all$ID = basename(IGV_all$ID)
write.table(IGV_all, file = 'Data/Final/IGV_all.seg', sep = '\t', row.names = F)
write.table(alterations, file = 'Data/Final/numericAlterations.txt', sep = '\t', row.names = F)




##----------------+
## transform numeric 
## alteration matrix in 
## binary format; no FILTERS considered
##----------------+
alterations = read.csv('Data/Final/numericAlterations.txt', sep = '\t')

## modify matrix
alteration_matrix = setNames(data.frame(matrix(ncol = 19, nrow = 0)), unique(alterations$gene))

all_out = data.frame()
counter = 1
for(id in unique(alterations$id)){
  data.mut.sub = alterations[which(alterations$id == id), ]
  if(nrow(data.mut.sub) != 0){
    for(j in unique(data.mut.sub$gene)){
      if(j %in% colnames(alteration_matrix)){
        alteration_matrix[counter, j] = data.mut.sub$n_call[which(data.mut.sub$gene == j)][1]
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

write.table(x = all_out, file = '~/Documents/MSKCC/Subhi/CSF/Data/Final/CSF_binary_noFilter.txt', sep = '\t')


##----------------+
## transform numeric
## alteration matrix in 
## binary format; consider
## only filter == PASS calls/genes
##----------------+
alterations = read.csv('Data/Final/numericAlterations.txt', sep = '\t')
alterations$n_call[which(alterations$filter %in% c('suppress_segment_too_large', 'suppress_large_homdel', 'suppress_likely_unfocal_large_gain'))] = NA

## modify matrix
alteration_matrix = setNames(data.frame(matrix(ncol = 19, nrow = 0)), unique(alterations$gene))

all_out = data.frame()
counter = 1
for(id in unique(alterations$id)){
  data.mut.sub = alterations[which(alterations$id == id), ]
  if(nrow(data.mut.sub) != 0){
    for(j in unique(data.mut.sub$gene)){
      if(j %in% colnames(alteration_matrix)){
        alteration_matrix[counter, j] = data.mut.sub$n_call[which(data.mut.sub$gene == j)][1]
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

write.table(x = all_out, file = '~/Documents/MSKCC/Subhi/CSF/Data/Final/CSF_binary_Filtered.txt', sep = '\t')


##----------------+
## only show FACETS QC TRUE
## samples in binary
##----------------+
sample_match = read.csv('Data/FINAL_samples/sample_match.txt', sep = '\t')
sample_match = sample_match[which(sample_match$fit == 'pass'), ]
sample_match = sample_match[!grepl(pattern = '.N0.*', x = sample_match$sample), ]

binary = read.csv('Data/Final/CSF_binary_Filtered.txt', sep = '\t')
colnames(binary) = gsub(pattern = '\\.', '-', colnames(binary))

binary_QC = binary[, which(colnames(binary) %in% sample_match$sample)]
write.table(binary_QC, file = 'Data/Final/CSF_binary_Filtered_QC_True.txt', sep = '\t', quote = F)



##-----------------
## work on ATRX:
## correct Sex assignment
## SampleSheet Excle!!!
##-----------------
# folders = list.files(path = '.')
# folders = folders[grepl(pattern = 'C-', x = folders)]
# sample_clinical = readxl::read_excel('Data/Final/ALEX_DATABASE_Nov6th_2022.xlsx')
# sample_clinical_short = sample_clinical[, c('Patient ID', 'Sample ID', 'Sex')]
# sample_clinical_short = as.data.frame(sample_clinical_short)
# missing_Sex = sample_clinical_short$`Patient ID`[which(sample_clinical_short$Sex == 'NA')]
Sex = read.csv('Data/FINAL_samples/Cohort_Sex.txt', sep = '\t')


# sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-K3KAAJ')] = 'F'
# sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-9DDAVH')] = 'M'
# sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-006880')] = 'F'
# sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-006881')] = 'M'
# sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-006882')] = 'M'
# sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-006883')] = 'M'
# sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-001521')] = 'M'
# sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-U2V48A')] = 'M'
# sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-AR2NUP')] = 'F'
# sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-002116')] = 'F'
# sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-0E4MEV')] = 'M'
# sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-001393')] = 'M'
# sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-WH70H5')] = 'M'
# sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-FXDRJ4')] = 'F'
# sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-K6YHWN')] = 'M'
# sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-6UMJWC')] = 'M'
# sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-F5JTRU')] = 'M'
# sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-07E77W')] = 'M'
# sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-J72UAD')] = 'F'
# sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-1K2R9U')] = 'F'
# sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-NWDR08')] = 'M'
# sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-9XEA6D')] = 'M'
# sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-ACM8V8')] = 'F'
# sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-JPXFN1')] = 'M'
# sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-WCM0NX')] = 'M'
# sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-4173R4')] = 'M'
# sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-23TLXA')] = 'M'
# sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-9E8T53')] = 'M'
# sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-N7LDHM')] = 'M'
# sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-002114')] = 'M'
# sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-V4TWV5')] = 'M'
# sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-R6KLRF')] = 'M'
# sample_clinical_short$Sex[which(sample_clinical_short$`Patient ID` == 'C-X4RMEV')] = 'F'

Sex$Sex[which(Sex$PatientID == 'C-001615')] = 'M'
Sex$Sex[which(Sex$PatientID == 'C-UMDEWL')] = 'M'
Sex$Sex[which(Sex$PatientID == 'C-DVCFA2')] = 'F'
Sex$Sex[which(Sex$PatientID == 'C-PTK9RU')] = 'M'
Sex$Sex[which(Sex$PatientID == 'C-D2P3RM')] = 'M'
Sex$Sex[which(Sex$PatientID == 'C-09T7CL')] = 'M'
Sex$Sex[which(Sex$PatientID == 'C-4173R4')] = 'F'
Sex$Sex[which(Sex$PatientID == 'C-56D79J')] = 'M'
Sex$Sex[which(Sex$PatientID == 'C-5XD5TE')] = 'M'
Sex$Sex[which(Sex$PatientID == 'C-65DNPK')] = 'M'
Sex$Sex[which(Sex$PatientID == 'C-6JMECR')] = 'M'
Sex$Sex[which(Sex$PatientID == 'C-9E8T53')] = 'M'
Sex$Sex[which(Sex$PatientID == 'C-9XPJL0')] = 'F'
Sex$Sex[which(Sex$PatientID == 'C-A32F9N')] = 'M'
Sex$Sex[which(Sex$PatientID == 'C-D7MUWF')] = 'M'
Sex$Sex[which(Sex$PatientID == 'C-NRMAAD')] = 'M'

# colnames(sample_clinical_short) = c('PatientID', 'SampleID', 'Sex')
# write.table(sample_clinical_short, file = 'Data/Final/Samples_Sex_Annotation.txt', sep = '\t', row.names = F)

#' out
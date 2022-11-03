##-----------------
## CCF clustering
##-----------------
##
## 10/31/2022
## chris-kreitzer

clean()
gc()
.rs.restartR()
setup(working.path = '~/Documents/MSKCC/Subhi/CSF/')


##-----------------
## Data
maf = read.csv('Data/FINAL_samples/maf_annotated.txt', sep = '\t')
csf_annotation = read.csv('Data/FINAL_samples/sample_match.txt', sep = '\t')
csf_clinical = readxl::read_excel('Data/FINAL_samples/CSF_Lastest_07102022.xlsx')
CSF = merge(maf, csf_annotation, by.x = 'Tumor_Sample_Barcode', by.y = 'sample', all.x = T)
CSF = merge(CSF, csf_clinical[,c('Sample ID', 'TYPE')], by.x = 'Tumor_Sample_Barcode', by.y = 'Sample ID', all.x = T)
binary = read.csv('Data/FINAL_samples/CSF_binary_all.txt', sep = '\t')
binary_true = read.csv('Data/FINAL_samples/CSF_binary_QC_true.txt', sep = '\t')

##-----------------
## matched pairs
##-----------------
sample_match = read.csv('Data/FINAL_samples/sample_match.txt', sep = '\t')
sample_original = readxl::read_excel('Data/FINAL_samples/CSF_Lastest_07102022.xlsx')
colnames(sample_original)[3] = 'SampleID'
colnames(sample_original)[2] = 'PatientID'

pass = merge(sample_match, sample_original, by.x = 'sample', by.y = 'SampleID', all.x = T)
pass = pass[!is.na(pass$TYPE) & !is.na(pass$ORDER), ]
pass = pass[which(pass$fit == 'pass'), ]


id_pairs = data.frame()
for(i in unique(pass$PATIENT_ID)){
  if(all(length(pass$sample[which(pass$PATIENT_ID == i)]) > 1,
     pass$TYPE[which(pass$PATIENT_ID == i & pass$ORDER == 1)] == 'TUMOR',
     pass$TYPE[which(pass$PATIENT_ID == i & pass$ORDER == 2)] == 'CSF')){
    id = i
    out = data.frame(id = id)
    
  } else next
  
  id_pairs = rbind(id_pairs, out)
}

sample_original = sample_original[which(sample_original$PatientID %in% id_pairs$id), ]


##----------------+
## solid tumor vs csf
##----------------+
gene_pairs = data.frame()
for(i in unique(id_pairs$id)){
  data_sub = sample_original[which(sample_original$PatientID == i), ]
  tumor_id = data_sub$SampleID[which(data_sub$ORDER == 1 & data_sub$TYPE == 'TUMOR')]
  tumor_id = ifelse(length(tumor_id) == 0, NA, tumor_id)
  csf_id = data_sub$SampleID[which(data_sub$ORDER == 2 & data_sub$TYPE == 'CSF')]
  csf_id = ifelse(length(csf_id) == 0, NA, csf_id)
  out = data.frame(id = i,
                   tumor = tumor_id,
                   csf = csf_id)
  gene_pairs = rbind(gene_pairs, out)
}

gene_pairs = gene_pairs[!with(gene_pairs, is.na(tumor) | is.na(csf)), ]


##----------------+
## CCF evolution
##----------------+
GOI = c('TP53', 'TERT', 'PTEN', 'IDH1', 'NF1', 'ATRX', 'EGFR', 'PIK3CA', 'RB1', 'CIC',
        'H3F3A', 'PDGFRA', 'BCOR', 'NOTCH2', 'PIK3R1', 'ARID1A', 'PTPN11', 'ATR', 'PMS2', 
        'FAT1', 'NOTCH1', 'TSC1', 'CDKN2A', 'CREBBP', 'FUBP1', 'BRAF')


a = data.frame(t = 0.4,
               c1 = 0.2,
               c2 = 0.5,
               c3 = 0.6)

ggplot(a, aes(x = c1, y = t)) +
  geom_point()
a



files = read.csv('~/Documents/MSKCC/Subhi/Lung/Facets_Chris.txt', sep = '\t')
head(files)
qc = read.csv('~/Documents/MSKCC/Subhi/Lung/Facets_All_qc.txt', sep = '\t')
head(qc)

all_out = data.frame()
for(i in 1:nrow(files)) {
  if(any(grepl(files$SAMPLE_ID[i], qc$qc))) {
    file = grep(pattern = files$SAMPLE_ID[i], x = qc$qc, value = T) 
  } else { 
    file = NA 
  } 
  
  out = data.frame(id = files$SAMPLE_ID[i], 
                   paths = file)
  
  all_out = rbind(all_out, out)
}

write.table(all_out, file = '~/Documents/MSKCC/Subhi/Lung/QC_paths_all.txt', sep = '\t', row.names = F)
write.table(QC_path, file = '~/Documents/MSKCC/Subhi/Lung/QC_paths_available.txt', sep = '\t', row.names = F)

datain = read.csv('~/Documents/MSKCC/Subhi/Lung/QC_all.txt', sep = '\t')
datain$fit = paste0(datain$path, datain$fit_name)

ids_keep = c()
for(i in unique(datain$tumor_sample_id)){
  if(length(datain$tumor_sample_id[which(datain$tumor_sample_id == i)]) > 1){
    data_sub = datain[which(datain$tumor_sample_id == i), ]
    if(any(data_sub$facets_qc)){
      fit = data_sub$fit[which(data_sub$facets_qc)]
    } else {
      fit = data_sub$fit[1]
    }
  } else {
    fit = datain$fit[which(datain$tumor_sample_id == i)]
  }
  ids_keep = c(ids_keep, fit)
}

QC_final = datain[which(datain$fit %in% ids_keep), ]
which(duplicated(QC_final$tumor_sample_id))

View(QC_final)


tumor_match = data.frame()
for(i in 1:nrow(gene_pairs)){
  tumor_maf = maf[which(maf$Tumor_Sample_Barcode == gene_pairs$tumor[i]), ]
  #csf_maf = maf[which(maf$Tumor_Sample_Barcode == gene_pairs$csf[i]), ]
  
  for(j in unique(GOI)){
    if(j %in% tumor_maf$Hugo_Symbol){
      gene = j
      type = 'Tumor'
      ccf = tumor_maf[which(tumor_maf$Hugo_Symbol == j), 'ccf_expected_copies']
      ccf = ifelse(length(ccf) > 1, ccf[1], ccf)
      clonality = tumor_maf[which(tumor_maf$Hugo_Symbol == j), 'clonality']
      clonality = ifelse(length(clonality) > 1, clonality[1], clonality)
      out = data.frame(id = gene_pairs$id[i],
                       gene = gene,
                       type = type,
                       ccf = ccf,
                       clonality = clonality)
    } else {
      gene = j
      type = 'Tumor'
      ccf = NA
      out = data.frame(id = gene_pairs$id[i],
                       gene = gene,
                       type = type,
                       ccf = ccf,
                       clonality = NA)
    } 
    
    # if(j %in% csf_maf$Hugo_Symbol){
    #   gene = j
    #   type = 'CSF'
    #   ccf = csf_maf[which(csf_maf$Hugo_Symbol == j), 'ccf_expected_copies']
    #   ccf = ifelse(length(ccf) > 1, ccf[1], ccf)
    #   out = data.frame(gene = gene,
    #                    type = type,
    #                    ccf = ccf)
    # } else {
    #   gene = j
    #   type = 'CSF'
    #   ccf = NA
    #   out = data.frame(gene = gene,
    #                    type = type,
    #                    ccf = ccf)
    # }
    tumor_match = rbind(tumor_match, out)
  }
}

##----------------+
## CSF match
csf_match = data.frame()
for(i in 1:nrow(gene_pairs)){
  #tumor_maf = maf[which(maf$Tumor_Sample_Barcode == gene_pairs$tumor[i]), ]
  csf_maf = maf[which(maf$Tumor_Sample_Barcode == gene_pairs$csf[i]), ]
  
  for(j in unique(GOI)){
    # if(j %in% tumor_maf$Hugo_Symbol){
    #   gene = j
    #   type = 'Tumor'
    #   ccf = tumor_maf[which(tumor_maf$Hugo_Symbol == j), 'ccf_expected_copies']
    #   ccf = ifelse(length(ccf) > 1, ccf[1], ccf)
    #   out = data.frame(id = gene_pairs$id[i],
    #                    gene = gene,
    #                    type = type,
    #                    ccf = ccf)
    # } else {
    #   gene = j
    #   type = 'Tumor'
    #   ccf = NA
    #   out = data.frame(id = gene_pairs$id[i],
    #                    gene = gene,
    #                    type = type,
    #                    ccf = ccf)
    # } 
    # 
    if(j %in% csf_maf$Hugo_Symbol){
      gene = j
      type = 'CSF'
      ccf = csf_maf[which(csf_maf$Hugo_Symbol == j), 'ccf_expected_copies']
      ccf = ifelse(length(ccf) > 1, ccf[1], ccf)
      clonality = csf_maf[which(csf_maf$Hugo_Symbol == j), 'clonality']
      clonality = ifelse(length(clonality) > 1, clonality[1], clonality)
      out = data.frame(id = gene_pairs$id[i],
                       gene = gene,
                       type = type,
                       ccf = ccf,
                       clonality = clonality)
    } else {
      gene = j
      type = 'CSF'
      ccf = NA
      out = data.frame(id = gene_pairs$id[i],
                       gene = gene,
                       type = type,
                       ccf = ccf,
                       clonality = NA)
    }
    
    csf_match = rbind(csf_match, out)
  }
}

tumor_csf_match = rbind(tumor_match, csf_match)


##' Visualization
ggplot(tumor_csf_match, aes(x = type, y = ccf, color = clonality)) +
  #geom_point() +
  geom_boxplot() +
  facet_grid(.~gene, space = 'free_x', scales = 'free_x')



  

csf_id = sample_original$SampleID[which(sample_original$PatientID == i & 
                                            sample_original$ORDER == 2)]
  
  tumor_maf = maf[which(maf$Tumor_Sample_Barcode == tumor_id), ]
  csf_maf = maf[which(maf$Tumor_Sample_Barcode == csf_id),]
  
  both = intersect(tumor_maf$Hugo_Symbol, csf_maf$Hugo_Symbol)
  
  if(length(both) != 0){
    tumor_ccf = tumor_maf[which(tumor_maf$Hugo_Symbol %in% both), 'ccf_expected_copies']
    csf_ccf = csf_maf[which(csf_maf$Hugo_Symbol %in% both), 'ccf_expected_copies']
    
    if(length(tumor_ccf) == length(csf_ccf) & length(tumor_ccf) == length(both)){
      out = data.frame(id = i,
                       gene = both,
                       tumor_ccf = tumor_maf[which(tumor_maf$Hugo_Symbol %in% both), 'ccf_expected_copies'],
                       csf_ccf = csf_maf[which(csf_maf$Hugo_Symbol %in% both), 'ccf_expected_copies'])
      
  } else next
  } else next

  b = rbind(b, out)
  # tumor_only = setdiff(tumor_maf$Hugo_Symbol, csf_maf$Hugo_Symbol)
  # csf_only = setdiff(csf_maf$Hugo_Symbol, tumor_maf$Hugo_Symbol)
}

ggplot(b, aes(x = tumor_ccf, y = csf_ccf)) +
  geom_jitter()



genes <- paste0("gene",1:1000)
set.seed(20210302)
gene_list <- list(A = sample(genes,100),
                  B = sample(genes,200),
                  C = sample(genes,300),
                  D = sample(genes,200))





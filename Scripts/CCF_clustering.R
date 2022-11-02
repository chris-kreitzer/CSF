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





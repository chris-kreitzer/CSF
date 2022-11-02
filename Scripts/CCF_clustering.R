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

maf = read.csv('Data/FINAL_samples/maf_annotated.txt', sep = '\t')
csf_annotation = read.csv('Data/FINAL_samples/sample_match.txt', sep = '\t')
csf_clinical = readxl::read_excel('Data/FINAL_samples/CSF_Lastest_07102022.xlsx')

CSF = merge(maf, csf_annotation, by.x = 'Tumor_Sample_Barcode', by.y = 'sample', all.x = T)
CSF = merge(CSF, csf_clinical[,c('Sample ID', 'TYPE')], by.x = 'Tumor_Sample_Barcode', by.y = 'Sample ID', all.x = T)
head(CSF)


boxplot(CSF$ccf_expected_copies[which(CSF$TYPE == 'TUMOR' & CSF$Hugo_Symbol == 'PTEN')],
     CSF$ccf_expected_copies[which(CSF$TYPE == 'CSF' & CSF$Hugo_Symbol == 'PTEN')])

head(csf_annotation$sample)
head(maf$Tumor_Sample_Barcode)



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

GOI = c("CDKN2A", "CDKN2B", "MTAP", 'EGFR', 'CDK4', 'PDGFRA', 'PTEN', 
        'KIT', 'MDM2', 'KDR', 'MDM4', 'RB1', 'MET', 'NF1', 'CDK6', 
        'TP53', 'KRAS', 'ATRX', 'FGF3', 'FGF4', 'FGF19')

a = data.frame()
for(i in unique(sample_original$PatientID)){
  if(length(sample_original$SampleID[which(sample_original$PatientID == i)]) > 1 &
     sample_original$TYPE[which(sample_original$ORDER == 1)] == 'TUMOR'){
    id = i
    out = data.frame(id = id)
  } else next
  a = rbind(a, out)
}

sample_original = sample_original[which(sample_original$PatientID %in% a$id), ]


##----------------+
## solid tumor vs csf
##----------------+
b = data.frame()
for(i in unique(a$id)){
  tumor_id = sample_original$SampleID[which(sample_original$PatientID == i & 
                                              sample_original$ORDER == 1)]
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

ggplot(b, aes(x = tumor_ccf[gene == 'TERT'], y = csf_ccf[gene == 'TERT'])) +
  geom_jitter()









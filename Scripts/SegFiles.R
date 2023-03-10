setwd('~/Documents/MSKCC/11_CSF/')
database = readxl::read_excel('00_Data/Database_Final+Feb24_ck_work_2023.xlsx', sheet = 1)

dim(database)
head(database)

samples_pass_tumor = database$`Sample ID`[which(database$Fit == 'full' & database$TYPE == 'TUMOR')]
samples_pass_csf = database$`Sample ID`[which(database$Fit == 'full' & database$TYPE == 'CSF')]

all_files = list.files(path = '01_countmatrices/', pattern = '.adjusted.seg', full.names = T, recursive = T)

seg_file_tumor = data.frame()
for(i in unique(samples_pass_tumor)){
  path = all_files[grep(pattern = i, x = all_files)]
  datain = read.csv(file = path, sep = '\t')
  datain$ID = i
  seg_file_tumor = rbind(seg_file_tumor, datain)
}

seg_file_csf = data.frame()
for(i in unique(samples_pass_csf)){
  path = all_files[grep(pattern = i, x = all_files)]
  datain = read.csv(file = path, sep = '\t')
  datain$ID = i
  seg_file_csf = rbind(seg_file_csf, datain)
}


write.table(x = seg_file_tumor, file = '00_Data/Segmentation_Tumor_QC_TRUE_03012023.txt', sep = '\t', row.names = F, quote = F)
write.table(x = seg_file_csf, file = '00_Data/Segmentation_CSF_QC_TRUE_03012023.txt', sep = '\t', row.names = F, quote = F)

cohort = readRDS('~/Documents/MSKCC/10_MasterThesis/Data/00_CohortData/Cohort_071322.rds')


##----------------+
## CCF annotate the gene calls
##----------------+

##-----------------
## CCF annotation of 
## mutations
##-----------------
##
## start: 11/02/2022
## revision: 11/11/2022
## revision: 03/06/2023
## revision: 06/18/2023
## 
## chris-kreitzer

clean()
gc()
.rs.restartR()
setup(working.path = '~/Documents/MSKCC/11_CSF/')
maf = read.csv('00_Data/MAF_oncokbanootated_June_18_2023.maf', sep = '\t')
database = readxl::read_excel('00_Data/CNA_Data_June_18_2023.xlsx', sheet = 1)
folders = list.files(path = '01_countmatrices/', full.names = T, recursive = F)


##-- MAF annotation only for QC FULL
full = database[which(database$Fit == 'full'), ]
retain = data.frame()
for(i in 1:nrow(full)){
  if(any(grepl(pattern = full$`Sample ID`[i], x = folders))){
    path = folders[grepl(pattern = full$`Sample ID`[i], x = folders)]
    out = data.frame(id = full$`Sample ID`[i],
                     path = path)
  } else next
  retain = rbind(retain, out)
}



##-----------------
## Annotation Function
##-----------------
all_out = data.frame()
for(i in 1:nrow(retain)){
  try({
    cncf = list.files(pattern = '_cncf.txt$', path = paste0(retain$path[i], '/'), full.names = T)
    cncf = read.csv(file = cncf, sep = '\t')
    qc = list.files(pattern = '.rds$', path = paste0(retain$path[i], '/'), full.names = T)
    qc = readRDS(file = qc)
    qc = facets_fit_qc(facets_output = qc)
    purity = qc$purity
    
    name = retain$id[i]
    print(name)
    
    if(name %in% maf$Tumor_Sample_Barcode){
      maf_sub = maf[which(maf$Tumor_Sample_Barcode == name), ]
      maf_annotated = facetsSuite::ccf_annotate_maf(maf = maf_sub,
                                                    segs = cncf,
                                                    purity = purity,
                                                    algorithm = 'em')
    } else next
    
    all_out = rbind(all_out, maf_annotated)
  })
}
    
 
write.table(all_out, file = '00_Data/MAF_annotated_03062023.txt', sep = '\t', row.names = F)



##-----------------
## MAF annotation NEW: June 18th, 2023
##-----------------

##-- DMP samples
dmp_full = database$id[which(database$CNA_fit == 'full')]
dmp_folders = list.files(path = '01_countmatrices/', full.names = T, recursive = F)

dmp_retain = data.frame()
for(i in 1:length(dmp_full)){
  if(any(grepl(pattern = dmp_full[i], x = dmp_folders))){
    path = dmp_folders[grepl(pattern = dmp_full[i], x = dmp_folders)]
    out = data.frame(id = dmp_full[i],
                     path = path)
  } else next
  dmp_retain = rbind(dmp_retain, out)
}
dmp_retain$Facets = 'single'


##-- Broad Run samples
facets_broad = database$id[grep(pattern = '*.Broad.*', x = database$CNA_fit)]
facets_folders = list.files(path = '07_CSF_refit/', full.names = T, recursive = F)

facets_retain = data.frame()
for(i in 1:length(facets_broad)){
  if(any(grepl(pattern = facets_broad[i], x = facets_folders))){
    path = facets_folders[grepl(pattern = facets_broad[i], x = facets_folders)]
    out = data.frame(id = facets_broad[i],
                     path = path)
  } else next
  facets_retain = rbind(facets_retain, out)
}
facets_retain$Facets = 'Broad'


##-- Hisens Run samples
facets_hisens = database$id[grep(pattern = '*.Hisens.*', x = database$CNA_fit)]
facets_folders = list.files(path = '07_CSF_refit/', full.names = T, recursive = F)

facets_retain_hi = data.frame()
for(i in 1:length(facets_hisens)){
  if(any(grepl(pattern = facets_hisens[i], x = facets_folders))){
    path = facets_folders[grepl(pattern = facets_hisens[i], x = facets_folders)]
    out = data.frame(id = facets_hisens[i],
                     path = path)
  } else next
  facets_retain_hi = rbind(facets_retain_hi, out)
}
facets_retain_hi$Facets = 'Hisens'


##-- MERGE
database_merged = rbind(dmp_retain, facets_retain, facets_retain_hi)


##-----------------
## Annotation Function
##-----------------
all_out = data.frame()
for(i in 1:nrow(database_merged)){
  try({
    #cncf = list.files(pattern = '_cncf.txt$', path = paste0(retain$path[i], '/'), full.names = T)
    #cncf = read.csv(file = cncf, sep = '\t')
    qc = list.files(pattern = '.rds$', path = paste0(database_merged$path[i], '/'), full.names = T)
    qc = ifelse(length(qc) > 1, qc[1], qc[1])
    qc = readRDS(file = qc)
    #qc = facets_fit_qc(facets_output = qc)
    #purity = qc$purity
    
    name = database_merged$id[i]
    print(name)
    
    if(name %in% maf$Tumor_Sample_Barcode){
      maf_sub = maf[which(maf$Tumor_Sample_Barcode == name), ]
      maf_annotated = facetsSuite::ccf_annotate_maf(maf = maf_sub,
                                                    segs = qc$segs,
                                                    purity = qc$purity,
                                                    algorithm = 'em')
    } else next
    
    all_out = rbind(all_out, maf_annotated)
  })
}

write.table(all_out, file = '00_Data/MAF_annotated_June_18_2023.txt', sep = '\t', row.names = F)


##----------------+
## Gene-Level Alteration all
##----------------+

# binary


#' out

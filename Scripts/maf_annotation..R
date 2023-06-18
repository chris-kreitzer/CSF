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
maf = read.csv('00_Data/Alex_Maf_oncokbanootated1_Feb_Hugo.maf', sep = '\t')
database = readxl::read_excel('00_Data/Database_Final+Feb24_ck_work_2023.xlsx', sheet = 1)
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


##----------------+
## Gene-Level Alteration all
##----------------+

# binary


#' out

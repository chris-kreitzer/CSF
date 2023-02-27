##----------------+
## create binary matrix
## for Heatmaps; numeric calls
##----------------+
##
## 02/25/2023
## chris-kreitzer


clean()
gc()
setwd('~/Documents/MSKCC/11_CSF/')

GOIs = c('CDKN2A','CDK4','CDK6','PTEN','EGFR','RET',
         'GATA3','PDGFRA','KIT','KDR','MET','MDM2',
         'MDM4','RB1','TP53')
database = readxl::read_excel('00_Data/Database_Final+Feb24_ck_work_2023.xlsx', sheet = 1)
database$id = basename(database$FacetsCountfile)

dmp = read.csv('00_Data/data_CNA.oncokb.txt.gz', sep = '\t')
dmp = dmp[,c('SAMPLE_ID', 'HUGO_SYMBOL', 'ALTERATION')]
dmp = dmp[which(dmp$HUGO_SYMBOL %in% GOIs), ]
dmp$ALTERATION = ifelse(dmp$ALTERATION == 'Deletion', '-2',
                        ifelse(dmp$ALTERATION == 'Amplification', '2', NA))
colnames(dmp) = c('id', 'gene', 'call')
dmp = dmp[which(dmp$id %in% database$`Sample ID`), ]
dmp$call = as.numeric(as.character(dmp$call))
full = read.csv('00_Data/FULL_FACETS_geneLevelSummary.txt', sep = '\t')
full = full[,c('id', 'gene', 'numeric_call')]
colnames(full)[3] = 'call'
reduced = read.csv('00_Data/REDUCED_FACETS_Sample_Summary.txt', sep = '\t')
reduced = reduced[,c('id', 'gene', 'call')]


facets_anno = rbind(full, reduced)
facets_anno_solo = facets_anno[!facets_anno$id %in% names(which(table(facets_anno$id) != 15)), ]
full_anno = data.frame()
for(i in names(which(table(facets_anno$id) != 15))){
  data_sub = facets_anno[which(facets_anno$id == i), ]
  genediff = setdiff(GOIs, data_sub$gene)
  data_sub = rbind(data_sub,
                   data.frame(id = i,
                              gene = genediff,
                              call = rep(0, length(genediff))))
  full_anno = rbind(full_anno, data_sub)
}

facets_annotated = rbind(full_anno, facets_anno_solo)
facets_annotated = merge(facets_annotated, database[,c('Sample ID', 'id')],
                         by.x = 'id', by.y = 'id', all.x = T)
colnames(facets_annotated)[4] = 'sample'


##----------------+
## final binary annotation
##----------------+
facets = facets_annotated
cbio = union(facets$sample[grep(pattern = '^P-0.*', x = facets$sample)], dmp$id)

csf = facets[!facets$sample %in% cbio, ]
csf = csf[,c('sample', 'gene', 'call')]
facetscbio = facets[which(facets$sample %in% cbio), ]

facetscbio = merge(facetscbio, dmp, by.x = c('sample', 'gene'), by.y = c('id','gene'), all.x = T)
colnames(facetscbio)[5] = 'cbio'
facetscbio$final = NA


facetscbio$final[which(facetscbio$call.x == facetscbio$cbio)] = facetscbio$cbio[which(facetscbio$call.x == facetscbio$cbio)]
facetscbio$final[which(is.na(facetscbio$call.x) & !is.na(facetscbio$cbio))] = facetscbio$cbio[which(is.na(facetscbio$call.x) & !is.na(facetscbio$cbio))]
facetscbio$final[which(!is.na(facetscbio$call.x) & is.na(facetscbio$cbio))] = facetscbio$call.x[which(!is.na(facetscbio$call.x) & is.na(facetscbio$cbio))]
facetscbio$final[which(facetscbio$call.x != facetscbio$cbio)] = facetscbio$cbio[which(facetscbio$call.x != facetscbio$cbio)]
facetscbio$final[which(is.na(facetscbio$call.x) & is.na(facetscbio$cbio))] = 0
facetscbio = facetscbio[,c('sample', 'gene', 'final')]


missing = data.frame()
for(i in unique(setdiff(cbio, unique(facetscbio$sample)))){
  out = data.frame(sample = i,
                   gene = GOIs,
                   call_x = NA)
  missing = rbind(missing, out)
}

missing = merge(missing, dmp, by.x = c('sample', 'gene'), by.y = c('id', 'gene'), all.x = T)
missing$final = NA
missing$final[which(is.na(missing$call_x) & !is.na(missing$call))] = missing$call[which(is.na(missing$call_x) & !is.na(missing$call))]
missing$final[which(is.na(missing$call_x) & is.na(missing$call))] = 0
missing = missing[,c('sample', 'gene', 'final')]

full_annotation = rbind(facetscbio, missing)
colnames(full_annotation)[3] = 'call'
cohort_out = rbind(full_annotation, csf)



##----------------+
## binary alteration matrix
##----------------+

## modify matrix
alteration_matrix = setNames(data.frame(matrix(ncol = 15, nrow = 0)), unique(cohort_out$gene))

all_out = data.frame()
counter = 1
for(id in unique(cohort_out$sample)){
  data.mut.sub = cohort_out[which(cohort_out$sample == id), ]
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

write.table(x = all_out, file = '~/Documents/MSKCC/11_CSF/00_Data/Binary_AlterationMatrix.txt', sep = '\t')



#' out
##-----------------
## CCF annotation of 
## mutations
##-----------------


clean()
gc()
.rs.restartR()
setup(working.path = '~/Documents/MSKCC/Subhi/CSF/')
maf = read.csv('Data/FINAL_samples/maf_onco_Alex.maf', sep = '\t')
maf = maf[which(maf$Tumor_Sample_Barcode != ""), ]

folders = list.files(path = '.')
folders = folders[grepl(pattern = 'C-', x = folders)]


##-----------------
## Annotation Function
##-----------------
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
        
        #' fetch values
        name = basename(j)
        print(name)
        segs = fit$cncf
        purity = fit$purity
        
        if(name %in% maf$Tumor_Sample_Barcode){
          maf_sub = maf[which(maf$Tumor_Sample_Barcode == name), ]
          maf_annotated = facetsSuite::ccf_annotate_maf(maf = maf_sub,
                                                        segs = segs,
                                                        purity = purity,
                                                        algorithm = 'em')
        } else next
        
      all_out = rbind(all_out, maf_annotated)
      }
    }
  })
}


write.table(all_out, file = 'Data/FINAL_samples/maf_annotated.txt', sep = '\t', row.names = F)

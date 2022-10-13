##---------------------------
## create binary Matrix for 
## CSF alterations:
##---------------------------
##
## start: 10/07/2022
## update: 10/12/2022
## chris-kreitzer
## 

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

alterations = data.frame()
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
        gene_out = gene_out[which(gene_out$gene %in% GOI), ]
        wgd = facets_fit_qc(facets_output = facets_out)$wgd
        name = basename(j)
        print(name)
        
        ##-----------
        ## gene-level-out:
        for(k in 1:nrow(gene_out)){
          if(wgd){
            copyStates = copynumberstates[which(copynumberstates$wgd == T), ]
            tcn = gene_out$tcn[k]
            lcn = gene_out$lcn[k]
            call = ifelse(is.na(lcn),
                          copyStates$numeric_call[which(copyStates$tcn == tcn &
                                                          is.na(copyStates$lcn))],
                          copyStates$numeric_call[which(copyStates$tcn == tcn &
                                                          copyStates$lcn == lcn)])
            call = ifelse(length(call) != 0, call, NA)
            gene = gene_out$gene[k]
            gene_final = data.frame(id = name,
                                    gene = gene,
                                    call = call)
          } else {
            copyStates = copynumberstates[which(copynumberstates$wgd == F), ]
            tcn = gene_out$tcn[k]
            lcn = gene_out$lcn[k]
            call = ifelse(is.na(lcn),
                          copyStates$numeric_call[which(copyStates$tcn == tcn &
                                                          is.na(copyStates$lcn))],
                          copyStates$numeric_call[which(copyStates$tcn == tcn &
                                                          copyStates$lcn == lcn)])
            call = ifelse(length(call) != 0, call, NA)
            gene = gene_out$gene[k]
            gene_final = data.frame(id = name,
                                    gene = gene,
                                    call = call)
          }
          alterations = rbind(alterations, gene_final)
        }
      }
      rm(fit, out, name, facets_out, gene_out, gene_final)
    }
  })
}

alterations = alterations[!is.na(alterations$gene), ]

##-----------------
## modify matrix
alteration_matrix = setNames(data.frame(matrix(ncol = 20,
                                               nrow = 0)), 
                                 unique(alterations$gene))

all_out = data.frame()
counter = 1

for(id in unique(alterations$id)){
  data.mut.sub = alterations[which(alterations$id == id), ]
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

write.table(x = all_out, file = '~/Documents/MSKCC/Subhi/CSF/CSF_test2.txt', sep = '\t')



##-----------------
## copy number states;
##-----------------

load('C-0E4MEV/s_C_0E4MEV_S011_d06/s_C_0E4MEV_S011_d06.Rdata')









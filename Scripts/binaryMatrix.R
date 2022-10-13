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

setup(working.path = '~/Documents/MSKCC/Subhi/CSF/')
folders = list.files(path = '.')
folders = folders[grepl(pattern = 'C-', x = folders)]

##' define alteration classes:
Amp = c("AMP (many states)", "AMP", "AMP (LOH)")
Gain = c("GAIN (many states)", "GAIN", "CNLOH & GAIN")
Loss = c("HETLOSS", "DOUBLE LOSS AFTER", "LOSS BEFORE", "CNLOH", 
         "CNLOH BEFORE & LOSS", "LOSS AFTER", "LOSS BEFORE & AFTER", "CNLOH AFTER", "LOSS & GAIN" )
Deletion = c("HOMDEL")
Diploid = c("DIPLOID", "DIPLOID or CNLOH", "TETRAPLOID")
Indeterminante = c('INDETERMINATE', NA)



copynumberstates = facetsSuite:::copy_number_states
GOI = c("CDKN2A", "CDKN2B", "MTAP", 'EGFR', 'CDK4', 'PDGFRA', 'PTEN', 
        'KIT', 'MDM2', 'KDR', 'MDM4', 'RB1', 'MET', 'NF1', 'CDK6', 
        'TP53', 'KRAS', 'ATRX', 'FGF3', 'FGF4', 'FGF19')

load(file = 'C-000499/P-0012463-T03-IM6/P-0012463-T03-IM6.Rdata')
fit$cncf$cf = NULL
fit$cncf$tcn = NULL
fit$cncf$lcn = NULL
fit$cncf = cbind(fit$cncf, cf = out$out$cf, tcn = out$out$tcn, lcn = out$out$lcn)
fit$cncf$lcn[fit$cncf$tcn == 1] = 0
fit$cncf$lcn.em[fit$cncf$tcn.em == 1] = 0

# Generate output
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

a = data.frame()
for(i in 1:nrow(gene_out)){
  if(wgd){
    copyStates = copynumberstates[which(copynumberstates$wgd == T), ]
    tcn = gene_out$tcn[i]
    lcn = gene_out$lcn[i]
    call = ifelse(is.na(lcn),
                  copyStates$numeric_call[which(copyStates$tcn == tcn &
                                                  is.na(copyStates$lcn))],
                  copyStates$numeric_call[which(copyStates$tcn == tcn &
                                                  copyStates$lcn == lcn)])
    call = ifelse(length(call) != 0, call, NA)
    gene = gene_out$gene[i]
    gene_final = data.frame(gene = gene,
                            call = call)
  } else {
    copyStates = copynumberstates[which(copynumberstates$wgd == F), ]
    tcn = gene_out$tcn[i]
    lcn = gene_out$lcn[i]
    call = ifelse(is.na(lcn),
                  copyStates$numeric_call[which(copyStates$tcn == tcn &
                                                  is.na(copyStates$lcn))],
                  copyStates$numeric_call[which(copyStates$tcn == tcn &
                                                  copyStates$lcn == lcn)])
    call = ifelse(length(call) != 0, call, NA)
    gene = gene_out$gene[i]
    gene_final = data.frame(gene = gene,
                            call = call)
  }
  a = rbind(a, gene_final)
}

a




all_out = data.frame()
for(i in unique(folders)){
  files = list.files(path = paste0(i, '/'), full.names = T, recursive = T)
  if(any(grepl(pattern = 'gene_level_out.txt', files))){
    path = files[grepl(pattern = 'gene_level_out.txt', files)]
    gene_level = read.csv(file = path, sep = '\t')
    
    #' create matrix
    alteration_matrix = setNames(data.frame(matrix(ncol = 14, 
                                                   nrow = 0)), 
                                 unique(gene_level$gene))
    
    
    counter = 1
    
    for(id in unique(gene_level$name)){
      data.mut.sub = gene_level[gene_level$name == id, ]
      if(nrow(data.mut.sub) != 0){
        for(j in unique(data.mut.sub$gene)){
          if(j %in% colnames(alteration_matrix)){
            alteration_matrix[counter, j] = data.mut.sub$copy_state[which(data.mut.sub$gene == j)][1]
            
          }
        }
      }
      
      alteration_matrix[counter, 'Sample.ID'] = id
      counter = counter + 1
      all_out = rbind(all_out, alteration_matrix)
    }
  } else next
  #all_out = rbind(all_out, alteration_matrix)
}

#' convert discrete to binary encoding
a = c(Amp, Gain, Loss, Deletion, Diploid)
a = as.character(unique(a))

for(i in unique(a)){
  if(i %in% Amp){
    all_out[all_out == i] = 2
  } else if (i %in% Gain) {
    all_out[all_out == i] = 1
  } else if (i %in% Diploid){
    all_out[all_out == i] = 0
  } else if (i %in% Loss) {
    all_out[all_out == i] = -1
  } else if (i %in% Deletion){
    all_out[all_out == i] = -2
  } else if (i %in% Indeterminante){
    all_out[all_out == i] = NA
  } else {
    next
  }
}

all_out[all_out == 'INDETERMINATE'] = NA
all_out = all_out[!is.na(all_out$CDKN2A), ]

all_out = t(all_out)
colnames(all_out) = all_out[nrow(all_out), ]
all_out = all_out[-nrow(all_out), ]

write.table(x = all_out, file = '~/Documents/MSKCC/Subhi/CSF/CSF_test.txt', sep = '\t')



##-----------------
## copy number states;
##-----------------

load('C-0E4MEV/s_C_0E4MEV_S011_d06/s_C_0E4MEV_S011_d06.Rdata')









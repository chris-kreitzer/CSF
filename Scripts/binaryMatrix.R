##---------------------------
## create binary Matrix for 
## CSF alterations:
##---------------------------
##
## 10/07/2022
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









## UTILITY functions for CSF analysis





## multi_snp-pileup read for FACETS analysis
## 
## start: 08/31/2022
## chris-kreitzer

clean()
gc()
# .rs.restartR()
setup(working.path = '~/Documents/MSKCC/Subhi/CSF/')


multi_readSnpMatrix = function(filename, 
                               tumor_sample = 2,
                               skip = 0L, 
                               err.thresh = 10, 
                               del.thresh = 10){
  
  pileup = read.csv(filename,
                    stringsAsFactors = FALSE)
  pileup[, 1] = as.character(pileup[,1])
  pileup[, 2] = as.numeric(as.integer(pileup[, 2]))
  pileup[, 3] = as.character(pileup[, 3])
  pileup[, 4] = as.character(pileup[, 4])
  cols = seq(5, length(pileup), 1)
  pileup[, cols] = apply(pileup[, cols], 2, function(x) as.numeric(as.integer(x)))
  
  if(tumor_sample - 1 > length(cols) - 4 / 4){
    stop(paste('wrong tumor number (column) selected'), call. = F)
  }
  
  #' select the right tumor sample:
  pileup_tumor = cbind(pileup[,c(1:8)], pileup[, grep(pattern = tumor_sample, x = colnames(pileup))])
  
  if(tumor_sample != 2){
    colnames(pileup_tumor)[9:12] = gsub(pattern = tumor_sample, replacement = '2', x = colnames(pileup_tumor)[9:12])
  }
  
  # remove chr if present in Chrom
  if(grepl("chr", pileup_tumor$Chromosome[1])){
    pileup_tumor$Chromosome = gsub("chr", "", pileup_tumor$Chromosome)
  }
  
  message(paste0((length(cols) - 4) / 4, ' tumor samples in pileup file; ready to analyze'))
  
  if(nrow(pileup_tumor) == 0){
    stop(paste(filename, 'does not exist or cannot be read properly.'), call. = F)
  }

  # remove loci where errors and deletions exceeded thresholds
  ii = which(pileup_tumor$File1E <= err.thresh & 
               pileup_tumor$File1D <= del.thresh & 
               pileup_tumor$File2E <= err.thresh & 
               pileup_tumor$File2D <= del.thresh)
  rcmat = pileup_tumor[ii, 1:2]
  rcmat$NOR.DP <- pileup_tumor$File1R[ii] + pileup_tumor$File1A[ii]
  rcmat$NOR.RD <- pileup_tumor$File1R[ii]
  rcmat$TUM.DP <- pileup_tumor$File2R[ii] + pileup_tumor$File2A[ii]
  rcmat$TUM.RD <- pileup_tumor$File2R[ii]
  
  return(rcmat)
}


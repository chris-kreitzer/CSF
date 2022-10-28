##-----------------
## chromosome X issue:
##-----------------

clean()
gc()
.rs.restartR()
source('~/Documents/GitHub/CSF/Scripts/UtilityFunctions.R')
setup(working.path = '~/Documents/MSKCC/Subhi/CSF/')
library(patchwork)

sample_sheet = readxl::read_excel('Data/FINAL_samples/CSF_Lastest_07102022.xlsx')
males = sample_sheet[which(sample_sheet$Sex == 'M'), c('Patient ID', 'Sample ID')]
females = sample_sheet[which(sample_sheet$Sex == 'F'), c('Patient ID', 'Sample ID')]

folders = list.files(path = '.')
folders = folders[grepl(pattern = 'C-', x = folders)]


##-----------------
## chromosome X issue:
##-----------------
chromosomeX = data.frame()
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
        name = basename(j)
        print(name)
        
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
        
        chromX = facets_out$segs[which(facets_out$segs$chrom == 23), ]
        chromX = paste(chromX$tcn.em, chromX$lcn.em, sep = ',', collapse = '|')
        
        out = data.frame(id = name,
                         cn_state = chromX)
        chromosomeX = rbind(chromosomeX, out)
      }
    }
  })
}

barplot(sort(table(chromosomeX$cn_state), decreasing = T),
        las = 2, cex.names = 0.9, cex.axis = 1, col = 'white')

par(mfrow = c(2,1))
#' Female
femaleX = chromosomeX[which(chromosomeX$id %in% females$`Sample ID`), ]
barplot(sort(table(femaleX$cn_state), decreasing = T),
        las = 2, cex.names = 0.9, cex.axis = 1, col = 'white')
mtext(text = 'female', side = 3, adj = 0.5)
#' male
maleX = chromosomeX[which(chromosomeX$id %in% males$`Sample ID`), ]
barplot(sort(table(maleX$cn_state), decreasing = T),
        las = 2, cex.names = 0.9, cex.axis = 1, col = 'white')
mtext(text = 'male', side = 3, adj = 0.5)


sample_match = read.csv('Data/FINAL_samples/sample_match.txt', sep = '\t')
sample_match = sample_match[which(sample_match$fit == 'pass'), ]


femaleX = femaleX[which(femaleX$id %in% sample_match$sample), ]


View(femaleX)
## C-002107
load('C-002116/P-0003822-T03-IM6/P-0003822-T03-IM6.Rdata')
fit$cncf$cf = NULL
fit$cncf$tcn = NULL
fit$cncf$lcn = NULL
fit$cncf = cbind(fit$cncf, cf = out$out$cf, tcn = out$out$tcn, lcn = out$out$lcn)
fit$cncf$lcn[fit$cncf$tcn == 1] = 0
fit$cncf$lcn.em[fit$cncf$tcn.em == 1] = 0

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

##----------------+
## Integer copy number
## correlation with CnLR
##----------------+

## start: 11/09/2022
## chris-kreitzer

clean()
gc()
.rs.restartR()
setup(working.path = '~/Documents/MSKCC/Subhi/CSF/')


## EGFR correlation: 
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

        #' fetch values
        name = basename(j)
        print(name)
        
        segs = fit$cncf
        purity = fit$purity
        
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
        
        EGFR_out = gene_out[which(gene_out$gene == 'EGFR'), c('median_cnlr_seg', 'tcn')]
        out = data.frame(id = name,
                         gene = 'EGFR',
                         tcn = EGFR_out$tcn,
                         cnlr = EGFR_out$median_cnlr_seg)
        
        all_out = rbind(all_out, out)
      }
    }
  })
}


## plot:
plot(all_out$tcn ~ all_out$cnlr,
     yaxt = 'n',
     ylab = '', 
     xlab = '')
abline(h = 6, lty = 'dashed', col = 'red')
axis(side = 2, at = seq(0, 180, 30), labels = seq(0, 180, 30), las = 2)
mtext(text = 'Integer copy number', side = 2, line = 2.5)
mtext(text = 'CopyNumber LogRatio', side = 1, line = 2.3)
mtext(text = 'EGFR', side = 3, line = 0.5, adj = 0, cex = 1.5)
text(x = -2.8, y = 160, labels = 'Pearson: 0.79\nSpearman: 0.83')

cor.test(all_out$tcn, all_out$cnlr, method = 'spearman')

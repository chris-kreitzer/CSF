DMP = c('P-0002647-T01-IM3','P-0002647-T02-IM5', 'P-0004837-T01-IM5','P-0005701-T01-IM5',
'P-0005852-T02-IM5','P-0007314-T01-IM5','P-0007835-T01-IM5','P-0009711-T01-IM5',
'P-0011213-T01-IM5','P-0013057-T01-IM5','P-0014099-T01-IM5','P-0016566-T01-IM6',
'P-0019261-T01-IM6','P-0019870-T01-IM6','P-0020652-T01-IM6','P-0021250-T01-IM6',
'P-0022938-T01-IM6','P-0023058-T01-IM6','P-0023638-T01-IM6','P-0028574-T01-IM6',
'P-0029256-T01-IM6','P-0030310-T01-IM6','P-0033342-T01-IM6','P-0034608-T01-IM6',
'P-0035158-T01-IM6','P-0036307-T01-IM6','P-0036307-T02-IM6','P-0036307-T04-IM6',
'P-0038169-T01-IM6','P-0041221-T01-IM6','P-0041221-T02-IM6','P-0042807-T01-IM6',
'P-0042807-T02-IM6','P-0042807-T03-IM6','P-0043097-T01-IM6','P-0043908-T01-IM6',
'P-0047557-T01-IM6','P-0048165-T01-IM6','P-0048165-T02-IM6','P-0048388-T01-IM6','P-0063096-T03-IM7','P-0063687-T01-IM7')

binary = read.csv('Data/FINAL_samples/CSF_binary_all.txt', sep = '\t')
colnames(binary) = gsub(pattern = '\\.', replacement = '-', x = colnames(binary))
clinical = readxl::read_excel('Data/FINAL_samples/CSF_Lastest_07102022.xlsx')

missing = data.frame()
for(i in 1:length(DMP)){
  if(DMP[i] %in% colnames(binary)){
    out = as.data.frame(binary[, DMP[i]], row.names = rownames(binary))
    colnames(out) = DMP[i]
  } else next
  missing = rbind(missing, out)
}

missing = binary[,which(colnames(binary) %in% DMP)]
not_missing = binary[,!colnames(binary) %in% DMP]
write.table(not_missing, file = '~/Desktop/notmissing.txt', sep = '\t')


load('C-10D913/P-0019261-T01-IM6/P-0019261-T01-IM6.Rdata')

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

gene = facetsSuite::gene_level_changes(facets_output = facets_out, genome = 'hg19', algorithm = 'em')
View(gene)
rm(facets_out, gene, fit, out)




binaryupdated = read.csv('~/Desktop/binaryMatrix_updated.txt', sep = '\t')
row.names(binaryupdated) = binaryupdated$gene
binaryupdated$gene = NULL
binaryupdated[binaryupdated == 1] = 2

write.table(binaryupdated, file = 'Data/FINAL_samples/CSF_binaryUPDATED.txt', sep = '\t')


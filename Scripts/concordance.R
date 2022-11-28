##-----------------
## concordance analysis:
##-----------------

clean()
gc()
.rs.restartR()
setup(working.path = '~/Documents/MSKCC/Subhi/CSF/')


##-----------------
sample_pairs = readxl::read_excel('Data/Final/SUBHI SPREADSHEET _USE.xlsx', sheet = 'DMP_CSF_Pairs')
sample_original = readxl::read_excel('Data/Final/SUBHI SPREADSHEET _USE.xlsx', sheet = 'Database')
sample_original_pass = sample_original[which(sample_original$FACETS_QC == 'pass'), ]
samples_alterations = read.csv('Data/Final/CSF_binary_Filtered_QC_True.txt', sep = '\t')
colnames(samples_alterations) = gsub(pattern = '\\.', replacement = '-', x = colnames(samples_alterations))
samples_alterations = samples_alterations[,-which(names(samples_alterations) %in% c('P-0006546-T02-IM5', 'P-0006546-T03-IM6',
                                                                                    's_C_006876_S013_d08', 's_C_HLCUF1_L001_d'))]

# sample_original = sample_original[, c('Sample ID', 'TYPE', 'ORDER')]
# colnames(sample_original)[1] = 'SampleID'
# 
# pass = merge(sample_match, sample_original, by.x = 'sample', by.y = 'SampleID', all.x = T)
# pass = pass[!is.na(pass$TYPE) & !is.na(pass$ORDER), ]
# pass = pass[which(pass$fit == 'pass'), ]
# 
# GOI = c("CDKN2A", "CDKN2B", "MTAP", 'EGFR', 'CDK4', 'PDGFRA', 'PTEN', 
#         'KIT', 'MDM2', 'KDR', 'MDM4', 'RB1', 'MET', 'NF1', 'CDK6', 
#         'TP53', 'KRAS', 'ATRX', 'FGF3', 'FGF4', 'FGF19')
# 
# alterations = read.csv('Data/FINAL_samples/CSF_cohort_purity_nalts.txt', sep = '\t')

alterations = data.frame()
for(i in 1:length(samples_alterations)){
  id = colnames(samples_alterations)[i]
  n_alts = sum(samples_alterations[,i] != 0, na.rm = T)
  n_deletions = sum(samples_alterations[,i] < 0, na.rm = T)
  n_amps = sum(samples_alterations[,i] > 0, na.rm = T)
  out = data.frame(id = id,
                   n_alts = n_alts,
                   n_deletions = n_deletions,
                   n_amps = n_amps)
  alterations = rbind(alterations, out)
}
write.table(alterations, file = 'Data/Final/n_Alterations.txt', sep = '\t', row.names = F)


##-----------------
## GENERAL overview: 
## Tumor vs CSFs
##-----------------
passed_samples = sample_original_pass[,c('Patient ID', 'Sample ID', 'TYPE', 'Purity', 'FACETS_QC', 'Oncogenic_CNA')]
passed_samples = as.data.frame(merge(passed_samples, alterations, by.x = 'Sample ID', by.y = 'id', all.x = T))



## number of all SCNA
dev.off()
#pdf(file = 'Figures/SCNA_all_comparison.pdf', width = 6, height = 6)
par(mfrow = c(1,3))
pdf(file = 'Figures/SCNA_FGA_Purity.pdf', width = 9, height = 6, paper = 'a4', onefile = T)

boxplot(passed_samples$n_alts[which(passed_samples$TYPE == 'TUMOR')],
        passed_samples$n_alts[which(passed_samples$TYPE == 'CSF')],
        xaxt = 'n',
        yaxt = 'n',
        ylim = c(0, 20))
axis(side = 1, at = c(1,2), labels = c('TUMOR', 'CSF'))
axis(side = 2, at = c(1, 5, 10 , 15, 20), labels = c(1, 5, 10 , 15, 20), las = 2)
box(lwd = 2)
mtext(text = paste0('p-value: ', round(t.test(passed_samples$n_alts[which(passed_samples$TYPE == 'TUMOR')],
                                              passed_samples$n_alts[which(passed_samples$TYPE == 'CSF')])$p.value, 3),
                    " (Welch's t-test)"), side = 3, line = 1.3)
mtext(text = '#SCNAs', side = 2, line = 2.8)
dev.off()



#' purity comparison
passed_samples$Purity = as.numeric(as.character(passed_samples$Purity))
pdf(file = 'Figures/purity_all_comparison.pdf', width = 6, height = 6)
boxplot(passed_samples$Purity[which(passed_samples$TYPE == 'TUMOR')],
        passed_samples$Purity[which(passed_samples$TYPE == 'CSF')],
        xaxt = 'n',
        yaxt = 'n',
        ylim = c(0, 1))
axis(side = 1, at = c(1,2), labels = c('TUMOR', 'CSF'))
axis(side = 2, at = c(0.1, 0.5, 1), labels = paste0(c(10, 50, 100), '%'), las = 2)
box(lwd = 2)
mtext(text = paste0('p-value: ', round(t.test(passed_samples$Purity[which(passed_samples$TYPE == 'TUMOR')],
                                              passed_samples$Purity[which(passed_samples$TYPE == 'CSF')])$p.value, 3),
                    " (Welch's t-test)"), side = 3, line = 1.3)
mtext(text = 'Purity', side = 2, line = 2.8)
dev.off()



##-----------------
## more metrics:
folders = list.files(path = '.')
folders = folders[grepl(pattern = 'C-', x = folders)]

all_out = data.frame()
for(i in unique(folders)){
  try({
    sub_dirs = list.dirs(path = i, full.names = T, recursive = F)
    if(length(sub_dirs) == 0) next
    else {
      for(j in unique(sub_dirs)){
        Rdata = list.files(pattern = '.Rdata$', path = paste0(j, '/'), full.names = T)
        load(file = Rdata)
        
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
        
        
        ##---- some metrics
        
        ploidy = facets_fit_qc(facets_output = facets_out)$ploidy
        wgd = facets_fit_qc(facets_output = facets_out)$wgd
        fga = facets_fit_qc(facets_output = facets_out)$fga
        n_amps = facets_fit_qc(facets_output = facets_out)$n_amps
        n_homdels = facets_fit_qc(facets_output = facets_out)$n_homdels
        n_loh = facets_fit_qc(facets_output = facets_out)$n_loh
        
        out = data.frame(name = name,
                           ploidy = ploidy,
                           wgd = wgd,
                           fga = fga,
                           n_amps = n_amps,
                           n_homdels = n_homdels,
                           n_loh = n_loh)

        all_out = rbind(all_out, out)
      }
    }
  })
}
write.table(all_out, file = 'Data/Final/FacetsSuite_n_Alterations.txt', sep = '\t', row.names = F)


alterations_all = merge(passed_samples, all_out, by.x = 'Sample ID', by.y = 'name', all.x = T)
head(alterations_all)


##-----------------
## FGA:
##-----------------
dev.off()
par(mfrow = c(1,3))
pdf(file = 'Figures/SCNA_FGA_Purity.pdf', width = 9, height = 6)
boxplot(alterations_all$fga[which(alterations_all$TYPE == 'TUMOR')],
        alterations_all$fga[which(alterations_all$TYPE == 'CSF')],
        xaxt = 'n',
        yaxt = 'n')
axis(side = 1, at = c(1, 2), labels = c('TUMOR', 'CSF'))
axis(side = 2, at = c(0.1, 0.5, 1), labels = paste0(c(10, 50, 100), '%'), las = 2)
box(lwd = 2)
mtext(text = paste0('p-value: ', round(t.test(alterations_all$fga[which(alterations_all$TYPE == 'TUMOR')],
                                              alterations_all$fga[which(alterations_all$TYPE == 'CSF')])$p.value, 3), 
                    " (Welch's t-test)"), side = 3, line = 1.3)
mtext(text = 'FGA', side = 2, line = 2.8)

dev.off()


##-----------------
## n_AMP
##-----------------
par(oma = c(3,3,3,3), omi = c(2,2,2,2))
pdf(file = 'Figures/nAMP_all_comparison.pdf', width = 6, height = 6)
boxplot(all_out$n_amps[which(all_out$type == 'TUMOR')],
        all_out$n_amps[which(all_out$type == 'CSF')],
        xaxt = 'n',
        yaxt = 'n',
        ylim = c(0, 10))
axis(side = 1, at = c(1, 2), labels = c('TUMOR', 'CSF'))
axis(side = 2, at = c(0, 5, 10), labels = c(1, 5, 10), las = 2)
box(lwd = 2)
mtext(text = paste0('p-value: ', round(t.test(all_out$n_amps[which(all_out$type == 'TUMOR')],
                                              all_out$n_amps[which(all_out$type == 'CSF')])$p.value, 3), ' (Welch Two Sample t-test)'), side = 3, line = 1.3 )
mtext(text = '#AMPs', side = 2, line = 2.8)

dev.off()

##-----------------
## n_homodeletions
##-----------------
par(oma = c(3,3,3,3), omi = c(2,2,2,2))
pdf(file = 'Figures/nHomoDel_all_comparison.pdf', width = 6, height = 6)
boxplot(all_out$n_homdels[which(all_out$type == 'TUMOR')],
        all_out$n_homdels[which(all_out$type == 'CSF')],
        xaxt = 'n',
        yaxt = 'n',
        ylim = c(0, 10))
axis(side = 1, at = c(1, 2), labels = c('TUMOR', 'CSF'))
axis(side = 2, at = c(0, 5, 10), labels = c(1, 5, 10), las = 2)
box(lwd = 2)
mtext(text = paste0('p-value: ', round(t.test(all_out$n_homdels[which(all_out$type == 'TUMOR')],
                                              all_out$n_homdels[which(all_out$type == 'CSF')])$p.value, 3), ' (Welch Two Sample t-test)'), side = 3, line = 1.3 )
mtext(text = '#HomoDels', side = 2, line = 2.8)

dev.off()


##-----------------
## n_LOH
##-----------------
par(oma = c(3,3,3,3), omi = c(2,2,2,2))
pdf(file = 'Figures/nLOH_all_comparison.pdf', width = 6, height = 6)
boxplot(all_out$n_loh[which(all_out$type == 'TUMOR')],
        all_out$n_loh[which(all_out$type == 'CSF')],
        xaxt = 'n',
        yaxt = 'n',
        ylim = c(0, 20))
axis(side = 1, at = c(1, 2), labels = c('TUMOR', 'CSF'))
axis(side = 2, at = c(0, 5, 10, 15, 20), labels = c(1, 5, 10, 15, 20), las = 2)
box(lwd = 2)
mtext(text = paste0('p-value: ', round(t.test(all_out$n_loh[which(all_out$type == 'TUMOR')],
                                              all_out$n_loh[which(all_out$type == 'CSF')])$p.value, 3), ' (Welch Two Sample t-test)'), side = 3, line = 1.3 )
mtext(text = '#LOH', side = 2, line = 2.8)

dev.off()


##-----------------
## concordance: 
## clonal alterations: 
## Tumor/CSFs PAIRS
##-----------------
sample_match = read.csv('Data/FINAL_samples/sample_match.txt', sep = '\t')
sample_original = readxl::read_excel('Data/FINAL_samples/CSF_Lastest_07102022.xlsx')
colnames(sample_original)[3] = 'SampleID'
colnames(sample_original)[2] = 'PatientID'

pass = merge(sample_match, sample_original, by.x = 'sample', by.y = 'SampleID', all.x = T)
pass = pass[!is.na(pass$TYPE) & !is.na(pass$ORDER), ]
pass = pass[which(pass$fit == 'pass'), ]

GOI = c("CDKN2A", "CDKN2B", "MTAP", 'EGFR', 'CDK4', 'PDGFRA', 'PTEN', 
        'KIT', 'MDM2', 'KDR', 'MDM4', 'RB1', 'MET', 'NF1', 'CDK6', 
        'TP53', 'KRAS', 'ATRX', 'FGF3', 'FGF4', 'FGF19')

a = data.frame()
for(i in unique(sample_original$PatientID)){
  if(length(sample_original$SampleID[which(sample_original$PatientID == i)]) > 1 &
     sample_original$TYPE[which(sample_original$ORDER == 1)] == 'TUMOR'){
    id = i
    out = data.frame(id = id)
  } else next
  a = rbind(a, out)
}

sample_original = sample_original[which(sample_original$PatientID %in% a$id), ]







##-----------------
## example: C-000499
## 1 solid tumor; 2 CSFs
##-----------------

load('C-000499/P-0012463-T03-IM6/P-0012463-T03-IM6.Rdata')
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


i = facetsSuite::cnlr_plot(facets_out, genome = 'hg19')
ii = facetsSuite::valor_plot(facets_out, genome = 'hg19')
iii = facetsSuite::icn_plot(facets_out, genome = 'hg19')
iv = facetsSuite::cf_plot(facets_out, genome = 'hg19')

i / ii / iii / iv

gene_level = facetsSuite::gene_level_changes(facets_output = facets_out, genome = 'hg19')
gene_level_goi = gene_level[which(gene_level$gene %in% GOI), ]

snps_solid = facets_out$snps

##-----------------
## CSF1:
##-----------------
load('C-000499/s_C_000499_L001_d/s_C_000499_L001_d.Rdata')
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


i = facetsSuite::cnlr_plot(facets_out, genome = 'hg19')
ii = facetsSuite::valor_plot(facets_out, genome = 'hg19')
iii = facetsSuite::icn_plot(facets_out, genome = 'hg19')
iv = facetsSuite::cf_plot(facets_out, genome = 'hg19')

i / ii / iii / iv

gene_level = facetsSuite::gene_level_changes(facets_output = facets_out, genome = 'hg19')
gene_level_goi = gene_level[which(gene_level$gene %in% GOI), ]

CSF1_snps = facets_out$snps


##-----------------
## CSF2:
##-----------------
load('C-000499/s_C_000499_L002_d/s_C_000499_L002_d.Rdata')
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


i = facetsSuite::cnlr_plot(facets_out, genome = 'hg19')
ii = facetsSuite::valor_plot(facets_out, genome = 'hg19')
iii = facetsSuite::icn_plot(facets_out, genome = 'hg19')
iv = facetsSuite::cf_plot(facets_out, genome = 'hg19')

i / ii / iii / iv

gene_level = facetsSuite::gene_level_changes(facets_output = facets_out, genome = 'hg19')
gene_level_goi = gene_level[which(gene_level$gene %in% GOI), ]
View(gene_level_goi)
View(gene_level)


##-----------------
## concordance clonal
## gene-wise
##-----------------

#' first rational:
#' whenever I see a clonal event in solid tumor tissue
#' I would expect to see the same alteration in CSF1,2,etc
GOI = c("CDKN2A", "CDKN2B", "MTAP", 'EGFR', 'CDK4', 'PDGFRA', 'PTEN', 
        'KIT', 'MDM2', 'KDR', 'MDM4', 'RB1', 'MET', 'NF1', 'CDK6', 
        'TP53', 'KRAS', 'FGF3', 'FGF4', 'FGF19')

all_out = data.frame()
for(dirs in list.dirs(path = 'C-000499/', recursive = F)){
  id = basename(dirs)
  load(file = list.files(pattern = '.Rdata', path = dirs, full.names = T))
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
  
  gene_level = facetsSuite::gene_level_changes(facets_output = facets_out, genome = 'hg19')
  gene_level_goi = gene_level[which(gene_level$gene %in% GOI), ]
  
  if(grepl(pattern = '^P-', x = id)){
    gene_level_clonal = gene_level_goi[which(gene_level_goi$filter %in% c('RESCUE', 'PASS') &
                                               gene_level_goi$cn_state != 'DIPLOID'), ]
    gene_level_clonal$clonality = ifelse(gene_level_clonal$cf.em >= facets_out$purity * 0.8, 'clonal', 'subclonal')
    gene_level_clonal = gene_level_clonal[which(gene_level_clonal$clonality == 'clonal'), ]
    n = length(unique(gene_level_clonal$gene))
    genes = paste(gene_level_clonal$gene, collapse = ',')
    
  } else {
    gene_level_clonal = gene_level_goi[which(gene_level_goi$filter %in% c('RESCUE', 'PASS') &
                                               gene_level_goi$cn_state != 'DIPLOID'), ]
    n = length(unique(gene_level_clonal$gene))
    genes = paste(gene_level_clonal$gene, collapse = ',')
  }
  
  out = data.frame(name = 'C-000499',
                   id = id,
                   genes = genes,
                   n = n,
                   order = NA)
  all_out = rbind(all_out, out)

}

all_out$concordance = NA
for(i in 1:nrow(all_out)){
  concordance = length(intersect(unlist(strsplit(all_out$genes[1], split = ',')), unlist(strsplit(all_out$genes[i], split = ','))))
  print(concordance)
  all_out$concordance[i] = concordance / all_out$n[1]
}


first = all_out
first$order = seq(1,3, 1)


##-----------------
all_out = data.frame()
for(dirs in list.dirs(path = 'C-000597//', recursive = F)){
  id = basename(dirs)
  load(file = list.files(pattern = '.Rdata', path = dirs, full.names = T))
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
  
  gene_level = facetsSuite::gene_level_changes(facets_output = facets_out, genome = 'hg19')
  gene_level_goi = gene_level[which(gene_level$gene %in% GOI), ]
  
  if(grepl(pattern = '^P-', x = id)){
    gene_level_clonal = gene_level_goi[which(gene_level_goi$filter %in% c('RESCUE', 'PASS') &
                                               gene_level_goi$cn_state != 'DIPLOID'), ]
    gene_level_clonal$clonality = ifelse(gene_level_clonal$cf.em >= facets_out$purity * 0.8, 'clonal', 'subclonal')
    gene_level_clonal = gene_level_clonal[which(gene_level_clonal$clonality == 'clonal'), ]
    n = length(unique(gene_level_clonal$gene))
    genes = paste(gene_level_clonal$gene, collapse = ',')
    
  } else {
    gene_level_clonal = gene_level_goi[which(gene_level_goi$filter %in% c('RESCUE', 'PASS') &
                                               gene_level_goi$cn_state != 'DIPLOID'), ]
    n = length(unique(gene_level_clonal$gene))
    genes = paste(gene_level_clonal$gene, collapse = ',')
  }
  
  out = data.frame(name = 'C-000597',
                   id = id,
                   genes = genes,
                   n = n,
                   order = NA)
  all_out = rbind(all_out, out)
  
}

all_out$concordance = NA
for(i in 1:nrow(all_out)){
  concordance = length(intersect(unlist(strsplit(all_out$genes[1], split = ',')), unlist(strsplit(all_out$genes[i], split = ','))))
  print(concordance)
  all_out$concordance[i] = concordance / all_out$n[1]
}

second = all_out
second$order = seq(1,3,1)


##-----------------
all_out = data.frame()
for(dirs in list.dirs(path = 'C-001393/', recursive = F)){
  id = basename(dirs)
  load(file = list.files(pattern = '.Rdata', path = dirs, full.names = T))
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
  
  gene_level = facetsSuite::gene_level_changes(facets_output = facets_out, genome = 'hg19')
  gene_level_goi = gene_level[which(gene_level$gene %in% GOI), ]
  
  if(grepl(pattern = '^P-', x = id)){
    gene_level_clonal = gene_level_goi[which(gene_level_goi$filter %in% c('RESCUE', 'PASS') &
                                               gene_level_goi$cn_state != 'DIPLOID'), ]
    gene_level_clonal$clonality = ifelse(gene_level_clonal$cf.em >= facets_out$purity * 0.8, 'clonal', 'subclonal')
    gene_level_clonal = gene_level_clonal[which(gene_level_clonal$clonality == 'clonal'), ]
    n = length(unique(gene_level_clonal$gene))
    genes = paste(gene_level_clonal$gene, collapse = ',')
    
  } else {
    gene_level_clonal = gene_level_goi[which(gene_level_goi$filter %in% c('RESCUE', 'PASS') &
                                               gene_level_goi$cn_state != 'DIPLOID'), ]
    n = length(unique(gene_level_clonal$gene))
    genes = paste(gene_level_clonal$gene, collapse = ',')
  }
  
  out = data.frame(name = 'C-001393/',
                   id = id,
                   genes = genes,
                   n = n,
                   order = NA)
  all_out = rbind(all_out, out)
  
}

all_out$concordance = NA
for(i in 1:nrow(all_out)){
  concordance = length(intersect(unlist(strsplit(all_out$genes[1], split = ',')), unlist(strsplit(all_out$genes[i], split = ','))))
  print(concordance)
  all_out$concordance[i] = concordance / all_out$n[1]
}

third = all_out
third$order = seq(1, 2, 1)


##-----------------
##-----------------
all_out = data.frame()
for(dirs in list.dirs(path = 'C-38YEPA', recursive = F)){
  id = basename(dirs)
  load(file = list.files(pattern = '.Rdata', path = dirs, full.names = T))
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
  
  gene_level = facetsSuite::gene_level_changes(facets_output = facets_out, genome = 'hg19')
  gene_level_goi = gene_level[which(gene_level$gene %in% GOI), ]
  
  if(grepl(pattern = '^P-', x = id)){
    gene_level_clonal = gene_level_goi[which(gene_level_goi$filter %in% c('RESCUE', 'PASS') &
                                               gene_level_goi$cn_state != 'DIPLOID'), ]
    gene_level_clonal$clonality = ifelse(gene_level_clonal$cf.em >= facets_out$purity * 0.8, 'clonal', 'subclonal')
    gene_level_clonal = gene_level_clonal[which(gene_level_clonal$clonality == 'clonal'), ]
    n = length(unique(gene_level_clonal$gene))
    genes = paste(gene_level_clonal$gene, collapse = ',')
    
  } else {
    gene_level_clonal = gene_level_goi[which(gene_level_goi$filter %in% c('RESCUE', 'PASS') &
                                               gene_level_goi$cn_state != 'DIPLOID'), ]
    n = length(unique(gene_level_clonal$gene))
    genes = paste(gene_level_clonal$gene, collapse = ',')
  }
  
  out = data.frame(name = 'C-38YEPA',
                   id = id,
                   genes = genes,
                   n = n,
                   order = NA)
  all_out = rbind(all_out, out)
  
}

all_out$concordance = NA
for(i in 1:nrow(all_out)){
  concordance = length(intersect(unlist(strsplit(all_out$genes[1], split = ',')), unlist(strsplit(all_out$genes[i], split = ','))))
  print(concordance)
  all_out$concordance[i] = concordance / all_out$n[1]
}

fifth = all_out
fifth$order = seq(1, 4, 1)



##-----------------
## Visualization:
##-----------------
all_out = rbind(first, second, third, fourth, fifth)
all_out$order = factor(all_out$order, levels = c(4,3,2,1))

ggplot(all_out, aes(x = name, y = order, fill = concordance)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) +
  scale_fill_gradient(
    low = "grey85",
    high = "red",
    space = "Lab",
    na.value = "grey50", name = 'CNA') +
  coord_fixed() +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  theme(panel.background = element_blank())


##-----------------
## all genes of interest:
##-----------------
Alterations_all = data.frame()
for(i in unique(pair_folders)){
  try({
    files = all_out$sample[which(all_out$PATIENT_ID == i)]
    sub_dirs = c(paste0('/Users/chriskreitzer/Documents/MSKCC/Subhi/CSF/', i, '/', files, '/'))
    #sub_dirs = list.dirs(path = i, full.names = T, recursive = F)
    
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
        
        gene_level = facetsSuite::gene_level_changes(facets_output = facets_out, genome = 'hg19')
        
        for(gene in GOI){
          print(gene)
          if(gene %in% unique(gene_level$gene)){
            gene_cnlr = gene_level[which(gene_level$gene == gene), 'median_cnlr_seg']
            gene_cnlr = ifelse(is.na(gene_cnlr), NA, gene_cnlr)
            out = data.frame(id = name,
                             gene = gene,
                             gene_cnlr = gene_cnlr)
          } else next
          Alterations_all = rbind(Alterations_all, out)
        }
      }
    }
  })
}


## modify output:
gene_list = list()
for(i in unique(Alterations_all$gene)){
  gene_sub = Alterations_all[which(Alterations_all$gene == i), ]
  gene_merge = merge(gene_sub, all_out, by.x = 'id', by.y = 'sample', all.x = T)
  gene_plot = data.frame()
  for(j in unique(gene_merge$PATIENT_ID)){
    tumor = gene_merge$gene_cnlr[which(gene_merge$PATIENT_ID == j & gene_merge$TYPE == 'TUMOR')]
    csf = gene_merge$gene_cnlr[which(gene_merge$PATIENT_ID == j & gene_merge$TYPE == 'CSF')]
    out = data.frame(id = j,
                     tumor = tumor,
                     csf = csf)
    gene_plot = rbind(gene_plot, out)
    
  }
  gene_list[[i]] = gene_plot
}

#' vis
plot_list = list()
for(i in seq_along(gene_list)){
  if(names(gene_list)[i] %in% c('CDKN2A', 'CDKN2B')){
    plot = ggplot(gene_list[[i]], aes(x = tumor, y = csf)) +
      geom_jitter() +
      geom_abline(slope = 1, intercept = 0, linetype = 'dashed') +
      scale_y_continuous(limits = c(-5, 1)) +
      scale_x_continuous(limits = c(-5, 1)) +
      theme(aspect.ratio = 1,
            panel.border = element_rect(fill = NA)) +
      labs(x = 'Tumor [CnLR]', y = '1.CSF [CnLR]', title = names(gene_list)[i]) +
      annotate(geom = 'text', x = -4, y = 0.7, 
               label = paste0('spearmans rho: ', 
                              round(cor.test(gene_list[[i]]$tumor, gene_list[[i]]$csf, method = 'spearman')$estimate[[1]], 3))) 
  } else {
    plot = ggplot(gene_list[[i]], aes(x = tumor, y = csf)) +
      geom_jitter() +
      geom_abline(slope = 1, intercept = 0, linetype = 'dashed') +
      scale_y_continuous(limits = c(-1, 4.5)) +
      scale_x_continuous(limits = c(-1, 4.5)) +
      theme(aspect.ratio = 1,
            panel.border = element_rect(fill = NA)) +
      labs(x = 'Tumor [CnLR]', y = '1.CSF [CnLR]', title = names(gene_list)[i]) +
      annotate(geom = 'text', x = 1.1, y = 4.1, 
               label = paste0('spearmans rho: ', 
                              round(cor.test(gene_list[[i]]$tumor, gene_list[[i]]$csf, method = 'spearman')$estimate[[1]], 3))) 
  }
  
  plot_list[[i]] = plot
  
}

library(cowplot)
plot_grid(plotlist = plot_list)


##----------------+
## IGV-like plot for
## paired CNA samples
##----------------+
sample_pairs = readxl::read_excel('Data/Final/SUBHI SPREADSHEET _USE.xlsx', sheet = 'DMP_CSF_Pairs')


##-- loop through folders
folders = list.files(path = '.')
folders = folders[grepl(pattern = 'C-', x = folders)]

IGV_all = data.frame()
for(i in unique(folders)){
  try({
    sub_dirs = list.dirs(path = i, full.names = T, recursive = F)
    if(length(sub_dirs) == 0) next
    else {
      for(j in unique(sub_dirs)){
        Rdata = list.files(pattern = '.Rdata$', path = paste0(j, '/'), full.names = T)
        load(file = Rdata)
        
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
        
        
        ##---- some metrics
        IGV = facetsSuite::format_igv_seg(facets_output = facets_out, normalize = T,
                                          sample_id = name)
        IGV_all = rbind(IGV_all, IGV)
      
      }
    }
  })
}

IGV_all$tag = NA

for(i in unique(IGV_all$ID)){
  if(i %in% sample_pairs$DMP_CNA_first){
    IGV_all$tag[which(IGV_all$ID == i)] = 'TUMOR'
  } else if (i %in% sample_pairs$CSF_CNA_first){
    IGV_all$tag[which(IGV_all$ID == i)] = 'CSF'
  } else {
    IGV_all$tag[which(IGV_all$ID == i)] = 'none'
  }
}

IGV_all = IGV_all[!IGV_all$tag %in% 'none', ]
write.table(IGV_all, file = 'Data/Final/IGV_paired_samples.txt', sep = '\t', row.names = F)





















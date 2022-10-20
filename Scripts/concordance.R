##-----------------
## concordance analysis:
##-----------------

clean()
gc()
.rs.restartR()
setup(working.path = '~/Documents/MSKCC/Subhi/CSF/')


##-----------------
sample_match = read.csv('Data/FINAL_samples/sample_match.txt', sep = '\t')
sample_original = readxl::read_excel('Data/FINAL_samples/CSF_Lastest_07102022.xlsx')
sample_original = sample_original[, c('Sample ID', 'TYPE', 'ORDER')]
colnames(sample_original)[1] = 'SampleID'

pass = merge(sample_match, sample_original, by.x = 'sample', by.y = 'SampleID', all.x = T)
pass = pass[!is.na(pass$TYPE) & !is.na(pass$ORDER), ]

GOI = c("CDKN2A", "CDKN2B", "MTAP", 'EGFR', 'CDK4', 'PDGFRA', 'PTEN', 
        'KIT', 'MDM2', 'KDR', 'MDM4', 'RB1', 'MET', 'NF1', 'CDK6', 
        'TP53', 'KRAS', 'ATRX', 'FGF3', 'FGF4', 'FGF19')




##-----------------
## GENERAL overview: Tumor vs CSFs






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


##-----------------
all_out = data.frame()
for(dirs in list.dirs(path = 'C-001121/', recursive = F)){
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
  
  out = data.frame(name = 'C-001121/',
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


##-----------------
all_out = data.frame()
for(dirs in list.dirs(path = 'C-001327/', recursive = F)){
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
  
  out = data.frame(name = 'C-001327/',
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




















all_out$order = seq(1, 3, 1)
test$order = seq(1, 4, 1)

ggplot(u, aes(x = name, y = rev(order), fill = n)) +
  geom_tile()

test = data.frame(name = 'C-0er',
                  id = c('a', 'b', 'c', 'd'),
                  n = c(9, 8 , 7, 3),
                  order = NA)

u = rbind(all_out, test)



gene_level_clonal = gene_level_goi[which(gene_level_goi$filter %in% c('RESCUE', 'PASS') &
                                           gene_level_goi$cn_state != 'DIPLOID'), ]
gene_level_clonal$clonality = ifelse(gene_level_clonal$cf.em >= fit$purity * 0.8, 'clonal', 'subclonal')
gene_level_clonal = gene_level_clonal[which(gene_level_clonal$clonality == 'clonal'), ]





con = merge(snps_solid, CSF1, by = 'maploc', all.x = T)
CDKN2A = con[which(con$chrom.x == 9 & con$maploc >= 21967751 & con$maploc <= 21995300), ]

plot(CDKN2A$cnlr.x, CDKN2A$cnlr.y)


head(con)
plot(con$cnlr.x, con$cnlr.y)
cor.test(CDKN2A$cnlr.x, CDKN2A$cnlr.y, method = 'spearman')

a = snps_solid[which(snps_solid$chrom == 9 & snps_solid$maploc >= 21967751 & snps_solid$maploc <= 21995300), ]
a$seq = seq(1, nrow(a), 1)
b = CSF1_snps[which(CSF1_snps$chrom == 9 & CSF1_snps$maploc >= 21967751 & CSF1_snps$maploc <= 21995300), ]
b$seq = seq(1, nrow(b), 1)

c = merge(a, b, by = 'seq', all.x = T)
plot(c$cnlr.x, c$cnlr.y)
cor.test(c$cnlr.x, c$cnlr.y)





##-----------------
## passed samples:
##-----------------
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



samples_pass = pass$sample[which(pass$fit %in% 'pass')]

purity_all = data.frame()
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
        
        purity = fit$purity
        out = data.frame(id = name,
                         purity = purity)
        purity_all = rbind(purity_all, out)
      }
    }
  })
}



n_alts_all = data.frame()
for(i in 1:length(alterations)){
  n = sum(alterations[, i] == 0, na.rm = T)
  n.na = sum(is.na(alterations[, i]))
  n_alts = 21 - sum(n, n.na)
  out = data.frame(id = colnames(alterations)[i],
                   n_alts = n_alts)
  n_alts_all = rbind(n_alts_all, out)
}

n_alts_all$id = gsub(pattern = '\\.', replacement = '-', n_alts_all$id)

str(purity_all)
str(n_alts_all)
pn = merge(purity_all, n_alts_all, by.x = 'id', by.y = 'id', all.x = T)
head(pn)
head(pn)
View(pn)

head(pn)
dim(pass)

pass = merge(pass, pn, by.x = 'sample', by.y = 'id', all.x = T)
dim(pass)
head(pass)
write.table(x = pass, file = '~/Desktop/CSF_cohort_purity_nalts.txt', sep = '\t', row.names = F)

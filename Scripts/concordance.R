##-----------------
## concordance analysis:
##-----------------

clean()
gc()
.rs.restartR()
setup(working.path = '~/Documents/MSKCC/Subhi/CSF/')


##----------------+
## just work on the 
## paired DMP/CSF samples
##----------------+
sample_pairs = readxl::read_excel('Data/Final/SUBHI SPREADSHEET _USE.xlsx', sheet = 'DMP_CSF_Pairs')
sample_original = readxl::read_excel('Data/Final/SUBHI SPREADSHEET _USE.xlsx', sheet = 'Database')
sample_original_pass = sample_original[which(sample_original$FACETS_QC == 'pass'), ]
samples_alterations = read.csv('Data/Final/CSF_binary_Filtered_QC_True.txt', sep = '\t')
colnames(samples_alterations) = gsub(pattern = '\\.', replacement = '-', x = colnames(samples_alterations))
samples_alterations = samples_alterations[,-which(names(samples_alterations) %in% c('P-0006546-T02-IM5', 'P-0006546-T03-IM6',
                                                                                    's_C_006876_S013_d08', 's_C_HLCUF1_L001_d'))]

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
## Tumor and matched CSF
##-----------------
sample_pairs = sample_pairs[, c('PatientID', 'DMP_CNA_first', 'CSF_CNA_first')]
passed_samples = sample_original_pass[,c('Patient ID', 'Sample ID', 'TYPE', 'Purity', 'FACETS_QC', 'Oncogenic_CNA')]
passed_samples = as.data.frame(merge(passed_samples, alterations, by.x = 'Sample ID', by.y = 'id', all.x = T))
passed_samples$plot = NA

for(i in unique(passed_samples$`Sample ID`)){
  if(i %in% sample_pairs$DMP_CNA_first){
    passed_samples$plot[which(passed_samples$`Sample ID` == i)] = 'DMP'
  } else if(i %in% sample_pairs$CSF_CNA_first){
    passed_samples$plot[which(passed_samples$`Sample ID` == i)] = 'CSF'
  } else {
    passed_samples$plot[which(passed_samples$`Sample ID` == i)] = 'none'
  }
}

passed_samples = passed_samples[!passed_samples$plot %in% 'none', ]


##----------------+
## FGA, n_amps and n_homo
##----------------+
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


##----------------+
## Plots; SCNA (n=20)
##----------------+
passed_samples = merge(passed_samples, all_out, by.x = 'Sample ID', by.y = 'name', all.x = T)
dev.off()
par(mfrow = c(1,3))
# pdf(file = 'Figures/SCNA_FGA_Purity.pdf', width = 9, height = 6, paper = 'a4', onefile = T)
boxplot(passed_samples$n_alts[which(passed_samples$plot == 'DMP')],
        passed_samples$n_alts[which(passed_samples$plot == 'CSF')],
        xaxt = 'n',
        yaxt = 'n',
        ylim = c(0, 20))
axis(side = 1, at = c(1,2), labels = c('TUMOR', 'CSF'))
axis(side = 2, at = c(1, 5, 10 , 15, 20), labels = c(1, 5, 10 , 15, 20), las = 2)
box(lwd = 2)
mtext(text = paste0('p-value: ', round(t.test(passed_samples$n_alts[which(passed_samples$plot == 'DMP')],
                                              passed_samples$n_alts[which(passed_samples$plot == 'CSF')])$p.value, 3),
                    " (Welch's t-test)"), side = 3, line = 1.3)
mtext(text = '#SCNAs', side = 2, line = 2.8)


##-- FGA
boxplot(passed_samples$fga[which(passed_samples$plot == 'DMP')],
        passed_samples$fga[which(passed_samples$plot == 'CSF')],
        xaxt = 'n',
        yaxt = 'n')
axis(side = 1, at = c(1, 2), labels = c('TUMOR', 'CSF'))
axis(side = 2, at = c(0.1, 0.5, 1), las = 2)
box(lwd = 2)
mtext(text = paste0('p-value: ', round(t.test(passed_samples$fga[which(passed_samples$plot == 'DMP')],
                                              passed_samples$fga[which(passed_samples$plot == 'CSF')])$p.value, 3), 
                    " (Welch's t-test)"), side = 3, line = 1.3)
mtext(text = 'Fraction Genome Altered', side = 2, line = 2.8)


##-- Purity
passed_samples$Purity = as.numeric(as.character(passed_samples$Purity))
boxplot(passed_samples$Purity[which(passed_samples$plot == 'DMP')],
        passed_samples$Purity[which(passed_samples$plot == 'CSF')],
        xaxt = 'n',
        yaxt = 'n',
        ylim = c(0, 1))
axis(side = 1, at = c(1,2), labels = c('TUMOR', 'CSF'))
axis(side = 2, at = c(0.1, 0.5, 1), labels = paste0(c(10, 50, 100), '%'), las = 2)
box(lwd = 2)
mtext(text = paste0('p-value: ', round(t.test(passed_samples$Purity[which(passed_samples$plot == 'DMP')],
                                              passed_samples$Purity[which(passed_samples$plot == 'CSF')])$p.value, 3),
                    " (Welch's t-test)"), side = 3, line = 1.3)
mtext(text = 'Purity', side = 2, line = 2.8)


##----------------+
## some other plots;
##----------------+
##' n_AMP
dev.off()
par(mfrow = c(1,3))
boxplot(passed_samples$n_amps.y[which(passed_samples$plot == 'DMP')],
        passed_samples$n_amps.y[which(passed_samples$plot == 'CSF')],
        xaxt = 'n',
        yaxt = 'n',
        ylim = c(0, 10))
axis(side = 1, at = c(1, 2), labels = c('TUMOR', 'CSF'))
axis(side = 2, at = c(0, 5, 10), labels = c(1, 5, 10), las = 2)
box(lwd = 2)
mtext(text = paste0('p-value: ', round(t.test(passed_samples$n_amps.y[which(passed_samples$plot == 'DMP')],
                                              passed_samples$n_amps.y[which(passed_samples$plot == 'CSF')])$p.value, 3), 
                    " (Welch's t-test)"), side = 3, line = 1.3)
mtext(text = '#AMPs; genome-wide', side = 2, line = 2.8)


##' n_homo-deletions
boxplot(passed_samples$n_homdels[which(passed_samples$plot == 'DMP')],
        passed_samples$n_homdels[which(passed_samples$plot == 'CSF')],
        xaxt = 'n',
        yaxt = 'n',
        ylim = c(0, 10))
axis(side = 1, at = c(1, 2), labels = c('TUMOR', 'CSF'))
axis(side = 2, at = c(0, 5, 10), labels = c(1, 5, 10), las = 2)
box(lwd = 2)
mtext(text = paste0('p-value: ', round(t.test(passed_samples$n_homdels[which(passed_samples$plot == 'DMP')],
                                              passed_samples$n_homdels[which(passed_samples$plot == 'CSF')])$p.value, 3), 
                    " (Welch's t-test)"), side = 3, line = 1.3)
mtext(text = '#Deep Deletions; genome-wide', side = 2, line = 2.8)


##' n_LOH
boxplot(passed_samples$n_loh[which(passed_samples$plot == 'DMP')],
        passed_samples$n_loh[which(passed_samples$plot == 'CSF')],
        xaxt = 'n',
        yaxt = 'n',
        ylim = c(0, 30))
axis(side = 1, at = c(1, 2), labels = c('TUMOR', 'CSF'))
axis(side = 2, at = c(0, 5, 10, 15, 20, 25, 30), labels = c(1, 5, 10, 15, 20, 25, 30), las = 2)
box(lwd = 2)
mtext(text = paste0('p-value: ', round(t.test(passed_samples$n_loh[which(passed_samples$plot == 'DMP')],
                                              passed_samples$n_loh[which(passed_samples$plot == 'CSF')])$p.value, 3), 
                    " (Welch's t-test)"), side = 3, line = 1.3)
mtext(text = '#LOH; genome-wide', side = 2, line = 2.8)


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























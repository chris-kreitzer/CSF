##----------------+
## How does the sequence coverage
## influence the gene-level calls in 
## CSF samples?
##----------------+

## start: 03/08/2023
## chris-kreitzer

clean()
setwd('~/Documents/MSKCC/11_CSF/')
source('~/Documents/GitHub/Y_chromosome_loss/PanCancer/Scripts/UtilityFunctions.R')
library(cowplot)
database = readxl::read_excel('00_Data/Database_Final_march12.xlsx')
CSF = list.files('04_MADSEQ/TumorCSF/', full.names = T)

##-- Read the Coverage Files
Coverage = data.frame()
for(i in unique(CSF)){
  print(i)
  datain = data.table::fread(input = i, sep = '\t')
  datain$ID = basename(i)
  Coverage = rbind(Coverage, datain)
}

Coverage$ID = substr(Coverage$ID, start = 1, stop = nchar(Coverage$ID) - 17)

ref_depth = data.frame()
for(i in unique(Coverage$ID)){
  ref = unique(Coverage$ref_depth[which(Coverage$ID == i)])
  id = i
  out = data.frame(id = id,
                   ref_depth = ref)
  ref_depth = rbind(ref_depth, out)
}

ref_depth = merge(ref_depth, database[,c('Sample.ID', 'TYPE', 'Gene.Panel', 'Fit')],
                  by.x = 'id', by.y = 'Sample.ID', all.x = T)

ref_depth = ref_depth[which(ref_depth$TYPE == 'CSF'), ]
ref_depth$Fit[which(ref_depth$Fit %in% c('N/Possible', 'N/Avail'))] = 'failed'
ref_depth$Fit[which(ref_depth$Fit %in% c('full', 'reduced'))] = 'passed'


##-- Investigate the output
CSF_plot = ggplot(ref_depth, 
       aes(x = ref_depth, color = Fit, fill = Fit)) +
  geom_histogram(bins = 100, color = 'black', linewidth = 0.2) +
  scale_fill_manual(values = c('failed' = 'red3',
                                'passed' = 'grey'), name = 'CNA-fitting') +
  theme_std(base_size = 14, base_line_size = 0.1) +
  scale_x_continuous(expand = c(0.01, 0),
                     limits = c(0, 2000)) +
  panel_border(size = 1, color = 'black') +
  scale_y_continuous(expand = c(0.0, 0.0),
                     limits = c(0, 40)) +
  labs(x = 'Average Seq. Depth', y = 'Number of cases', title = 'CSF-Tumor Samples')

CSF_plot

ggsave_golden(filename = '05_Plots/RefDepth_CSF.pdf', plot = CSF_plot, width = 8)


##-- TUMOR samples:
tumor = read.csv('00_Data/Sample_coverage.txt', sep = '\t')
tumor = merge(tumor, database[,c('Sample.ID', 'TYPE', 'Gene.Panel', 'Fit')],
              by.x = 'Sample.ID', by.y = 'Sample.ID', all.x = T)

tumor = tumor[which(tumor$TYPE == 'TUMOR'), ]
tumor$Fit[which(tumor$Fit %in% c('N/Possible', 'N/Avail', 'cbio'))] = 'failed'
tumor$Fit[which(tumor$Fit %in% c('full', 'reduced'))] = 'passed'


tumor_plot = ggplot(tumor, 
                  aes(x = Sample.coverage, color = Fit, fill = Fit)) +
  geom_histogram(bins = 100, color = 'black', linewidth = 0.2) +
  scale_fill_manual(values = c('failed' = 'red3',
                               'passed' = 'grey'), name = 'CNA-fitting') +
  theme_std(base_size = 14, base_line_size = 0.1) +
  scale_x_continuous(expand = c(0.01, 0),
                     limits = c(0, 1500)) +
  panel_border(size = 1, color = 'black') +
  scale_y_continuous(expand = c(0.0, 0.0),
                     limits = c(0, 12)) +
  labs(x = 'Average Seq. Depth', y = 'Number of cases', title = 'Solid Tumor Samples')

tumor_plot
ggsave_golden(filename = '05_Plots/RefDepth_Tumor.pdf', plot = tumor_plot, width = 8)










##-- how much of the samples with average seq. depth <500 
## pass QC?

out = data.frame(SAMPLE_ID = 's_C_001327_S033_d10_normal',
                 GENE_PANEL = 'IMPACT468',
                 Normal_bam = '/juno/res/dmpcollab/dmpshare/share/irb12_245/Y/F/YF000268-N.bam')

write.table(x = out, file = '~/Desktop/normal.txt', sep = '\t', row.names = F, quote = F)


cou = read.csv('01_countmatrices/s_C_001327_N901_dZ_IM6.rg.md.abra.printreads__s_C_001327_S033_d10.rg.md.abra.printreads.dat.gz/s_C_001327_N901_dZ_IM6.rg.md.abra.printreads__s_C_001327_S033_d10.rg.md.abra.printreads.dat.gz', sep = ',')
t.test(cou$File1R, cou$File2R, paired = T)


aa = Coverage[which(Coverage$ref_depth < 50), ]



axis(side = 2, las = 2, line = -1)
title(main = 's_C_000499_L001', adj = 0)
title(main = 'CSF-Tumor sample; normed depth', adj = 1)
dev.off()



b = Coverage[which(Coverage$ID == 's_C_000624_L001_d_normed_depth.txt'), ]
plot(hist(b$normed_depth, nclass = 50),
     yaxt = 'n',
     xlab = '',
     main = '')
abline(v = unique(b$ref_depth), col = 'red', lwd = 2)
axis(side = 2, las = 2, line = -1)



plo = read.csv('00_Data/Database_Final+Feb24_ck_work_2023.txt', sep = '\t')
plo = plo[,c('Sample.ID', 'ploidy')]
plo$ploidy[which(plo$ploidy == 'N.A')] = NA
plo$ploidy = as.numeric(as.character(plo$ploidy))


write.table(x = plo, file = '~/Desktop/CSF_ploidies.txt', sep = '\t', row.names = F, quote = F)

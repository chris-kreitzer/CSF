##-----------------
## Genome-wide plot (CSF VS TUMOR)
## 07/07/2023
##-----------------
library(cowplot)
library(ggplot2)
library(tidyverse)

database = readxl::read_excel('00_Data/CNA_Data_June_18_2023.xlsx')
databaseTRUE = database[which(database$QC == TRUE), ]


csf = list.files(path = '07_CSF_refit/', pattern = '*.pass.rds$', full.names = T, recursive = T)
tumor = list.files(path = '01_countmatrices/', pattern = '*.rds$', full.names = T, recursive = T)

databaseCSF = databaseTRUE[which(databaseTRUE$TYPE == 'CSF'), ]
databaseTUMOR = databaseTRUE[which(databaseTRUE$TYPE == 'Tumor'), ]


CSF_IGV = data.frame()
for(i in 1:nrow(databaseCSF)){
  if(length(grep(pattern = databaseCSF$id[i], x = csf)) == 1){
    datain = readRDS(file = csf[grep(pattern = databaseCSF$id[i], x = csf)])
    dataconverted = facetsSuite::format_igv_seg(facets_output = datain, sample_id = basename(databaseCSF$id[i]), normalize = T)
  }
  else {
    files = csf[grep(pattern = databaseCSF$id[i], x = csf)]
    file_use = files[grep(pattern = '*._third_pass.rds$', x = files)]
    datain = readRDS(file = file_use)
    dataconverted = facetsSuite::format_igv_seg(facets_output = datain, sample_id = basename(databaseCSF$id[i]), normalize = T)
  }
  CSF_IGV = rbind(CSF_IGV, dataconverted)
} 


##-- Tumor
TUMOR_IGV = data.frame()

for(i in 1:nrow(databaseTUMOR)){
  if(length(grep(pattern = databaseTUMOR$id[i], x = tumor)) == 1){
    datain = readRDS(file = tumor[grep(pattern = databaseTUMOR$id[i], x = tumor)])
    dataconverted = facetsSuite::format_igv_seg(facets_output = datain, sample_id = basename(databaseTUMOR$id[i]), normalize = T)
  }
  else {print(databaseTUMOR$id[i])}
  TUMOR_IGV = rbind(TUMOR_IGV, dataconverted)
} 




##-----------------
## Plot function

## convert CSF data
CSF.Data = CSF_IGV
colnames(CSF.Data)[6] = 'LRR'
CSF.Data$group = ifelse(CSF.Data$LRR < -0.2, 'loss',
                           ifelse(CSF.Data$LRR > -0.2 & CSF.Data$LRR < 0.2, 'neutral', 'gain'))

## convert Tumor data
Tumor.Data = TUMOR_IGV
Tumor.Data$group = ifelse(Tumor.Data$seg.mean < -0.2, 'loss',
                         ifelse(Tumor.Data$seg.mean > -0.2 & Tumor.Data$seg.mean < 0.2, 'neutral', 'gain'))


###################
## binning function:
bin.function = function(data, bin.size = 100000000){
  Facets = data
  Facets = Facets[Facets$chrom %in% c(1:23),, drop = F]
  
  bins = data.frame()
  for(i in unique(Facets$chrom)){
    chrom.min = min(Facets$loc.start[which(Facets$chrom == i)]) 
    chrom.max = max(Facets$loc.end[which(Facets$chrom == i)]) 
    bin.seq = seq(chrom.min, chrom.max, bin.size)
    
    roi = data.frame(bin = paste0(seq(1, length(bin.seq), by = 1), 'bin'),
                     start = bin.seq, 
                     end = c(bin.seq[-1], chrom.max))
    
    roi$start = roi$start + 1
    roi$chrom = i
    
    bins = rbind(bins, roi)
  }
  
  return(bins)
}


###################
## segmentation algorithm

segmentation.function = function(data, group, bin.size, chromosomes = c(1:23)){
  message('loss; neutral and gain groups are available to select')
  sample.size = length(unique(data$ID))
  # select dataset (either Facets, or TCGA)
  # selcet group: gain or loss
  Facets.selected = data[which(data$group == as.character(group)),, drop = F]
  Facets.selected = Facets.selected[Facets.selected$chrom %in% chromosomes,, drop = F]
  
  # run bin.function
  bins = bin.function(data = data, bin.size = bin.size)
  
  out = data.frame()
  for(i in unique(bins$chrom)){
    bins.sub = bins[bins$chrom == i, ]
    bins.sub = bins.sub[order(bins.sub$start), ]
    facets.sub = Facets.selected[Facets.selected$chrom == i, ]
    
    for(alteration in 1:nrow(bins.sub)){
      segment = facets.sub[which(as.numeric(facets.sub$loc.start) <= as.numeric(bins.sub[alteration, 'start']) &
                                   as.numeric(facets.sub$loc.end) >= as.numeric(bins.sub[alteration, 'end']) |
                                   as.numeric(facets.sub$loc.start) >= as.numeric(bins.sub[alteration, 'start']) &
                                   as.numeric(facets.sub$loc.start) < as.numeric(bins.sub[alteration, 'end']) &
                                   as.numeric(facets.sub$loc.end) >= as.numeric(bins.sub[alteration, 'end']) |
                                   as.numeric(facets.sub$loc.start) <= as.numeric(bins.sub[alteration, 'start']) &
                                   as.numeric(facets.sub$loc.end) <= as.numeric(bins.sub[alteration, 'end']) &
                                   as.numeric(facets.sub$loc.end) > as.numeric(bins.sub[alteration, 'start'])), ]
      
      # print(length(unique(segment$ID)))
      ratio = length(unique(segment$ID)) / sample.size 
      group = segment$group
      
      x = bins.sub[alteration, ]
      x$ratio = ratio
      x$chrom = i
      
      out = rbind(out, x)
      
    }
    
    rm(i, bins.sub, facets.sub, segment, ratio, group, x)
    
  }
  
  return(out)
  
}


#######################################
## concentrate on gained regions first
csf.gain = segmentation.function(data = CSF.Data, group = 'gain', bin.size = 1000000)
csf.gain = csf.gain[!csf.gain$bin %in% c('248bin', '249bin', '250bin'), ]
csf.gain = csf.gain[!(csf.gain$chrom == 7 & csf.gain$bin == '160bin'), ]
csf.gain = csf.gain[!(csf.gain$chrom == 9 & csf.gain$bin == '142bin'), ]
csf.gain = csf.gain[!(csf.gain$chrom == 12 & csf.gain$bin %in% c('134bin', '133bin')), ]
csf.gain = csf.gain[!(csf.gain$chrom == 14 & csf.gain$bin %in% c('89bin', '88bin')), ]
csf.gain = csf.gain[!(csf.gain$chrom == 15 & csf.gain$bin %in% c('83bin', '82bin', '81bin')), ]
csf.gain = csf.gain[!(csf.gain$chrom == 16 & csf.gain$bin %in% c('91bin', '90bin')), ]
csf.gain = csf.gain[!(csf.gain$chrom == 17 & csf.gain$bin %in% c('82bin', '81bin', '80bin')), ]
csf.gain = csf.gain[!(csf.gain$chrom == 18 & csf.gain$bin %in% c('78bin', '77bin')), ]
csf.gain = csf.gain[!(csf.gain$chrom == 20 & csf.gain$bin %in% c('63bin')), ]
csf.gain = csf.gain[!(csf.gain$chrom == 21 & csf.gain$bin %in% c('39bin')), ]
csf.gain = csf.gain[!(csf.gain$chrom == 22 & csf.gain$bin %in% c('35bin')), ]
csf.gain = csf.gain[!(csf.gain$chrom == 23), ]

tumor.gain = segmentation.function(data = Tumor.Data, group = 'gain', bin.size = 1000000)
tumor.gain = tumor.gain[!(tumor.gain$chrom == 23), ]




combined.gain.data = cbind(csf.gain, tumor.gain[, c('ratio')])
colnames(combined.gain.data)[5] = 'ratio.csf'
colnames(combined.gain.data)[6] = 'ratio.tumor'
combined.gain.data$bin = paste0(seq(1, nrow(combined.gain.data), 1), 'bin')
combined.gain.data$bin = factor(combined.gain.data$bin, levels = combined.gain.data$bin)

## fetch chromosome ends
chrom.breaks = combined.gain.data %>% arrange(chrom, end) %>% group_by(chrom) %>% summarise(max(end))
chrom.breaks$chrom = as.integer(as.character(chrom.breaks$chrom))
chrom.breaks = chrom.breaks[order(chrom.breaks$chrom), ]
colnames(chrom.breaks)[2] = 'end'

chrom.out = combined.gain.data$bin[which(combined.gain.data$end %in% chrom.breaks$end)]


## add chromosome labels on top of gain plot
x = extract_numeric(chrom.out)
x = c(0, x)
diff.all = c()
for(i in 1:length(x)){
  diff = (x[i+1] - x[i])
  diff.all = c(diff.all, diff)
}
diff.all = diff.all[-length(diff.all)]
diff.quote = diff.all / 2
x = x[-1]
diff.quote

bin.ve = paste0(ceiling(x - diff.quote), 'bin')
bin.ve = factor(bin.ve, levels = bin.ve)


#######################################
## make plot
## select the two color bars; either way 
csf.color = '#8073AC'
tumor.color = '#D6604D'
cancer.type = 'ctDNA in CSF (Glioma)'
sample.size = length(unique(databaseTRUE$id))


gain.plot = ggplot(combined.gain.data, aes(x = bin)) +
  geom_bar(aes(y = ratio.tumor), stat = 'identity', alpha = 0.5, width = 1, fill = tumor.color) +
  geom_bar(aes(y = ratio.csf), stat = 'identity', alpha = 1, width = 1, fill = csf.color) +
  
  geom_vline(xintercept = chrom.out, linetype = "dashed", colour = "grey") +
  geom_path(aes(y = ratio.tumor), group = 1, col = 'black') +
  theme_cowplot() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        panel.background = element_blank(),
        legend.position = "right",
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  scale_x_discrete(drop = FALSE) +
  ylim(0, 1) +
  ylab("Prop. of Patients\nWith Gain/Amp") +
  
  annotate('text',
           x = bin.ve,
           y = 1,
           label = seq(1, length(bin.ve), 1)) +
  labs(title = paste0(cancer.type, ' ', '(n=', sample.size, ')'))




##-- LOSS plots
csf.loss = segmentation.function(data = CSF.Data, group = 'loss', bin.size = 1000000)
csf.loss = csf.loss[!csf.loss$bin %in% c('248bin', '249bin', '250bin'), ]
csf.loss = csf.loss[!(csf.loss$chrom == 7 & csf.loss$bin == '160bin'), ]
csf.loss = csf.loss[!(csf.loss$chrom == 9 & csf.loss$bin == '142bin'), ]
csf.loss = csf.loss[!(csf.loss$chrom == 12 & csf.loss$bin %in% c('134bin', '133bin')), ]
csf.loss = csf.loss[!(csf.loss$chrom == 14 & csf.loss$bin %in% c('89bin', '88bin')), ]
csf.loss = csf.loss[!(csf.loss$chrom == 15 & csf.loss$bin %in% c('83bin', '82bin', '81bin')), ]
csf.loss = csf.loss[!(csf.loss$chrom == 16 & csf.loss$bin %in% c('91bin', '90bin')), ]
csf.loss = csf.loss[!(csf.loss$chrom == 17 & csf.loss$bin %in% c('82bin', '81bin', '80bin')), ]
csf.loss = csf.loss[!(csf.loss$chrom == 18 & csf.loss$bin %in% c('78bin', '77bin')), ]
csf.loss = csf.loss[!(csf.loss$chrom == 20 & csf.loss$bin %in% c('63bin')), ]
csf.loss = csf.loss[!(csf.loss$chrom == 21 & csf.loss$bin %in% c('39bin')), ]
csf.loss = csf.loss[!(csf.loss$chrom == 22 & csf.loss$bin %in% c('35bin')), ]
csf.loss = csf.loss[!(csf.loss$chrom == 23), ]

tumor.loss = segmentation.function(data = Tumor.Data, group = 'loss', bin.size = 1000000)
tumor.loss = tumor.loss[!(tumor.loss$chrom == 23), ]


combined.loss.data = cbind(csf.loss, tumor.loss[, c('ratio')])
colnames(combined.loss.data)[5] = 'ratio.csf'
colnames(combined.loss.data)[6] = 'ratio.tumor'
combined.loss.data$bin = paste0(seq(1, nrow(combined.loss.data), 1), 'bin')
combined.loss.data$bin = factor(combined.loss.data$bin, levels = combined.loss.data$bin)

## fetch chromosome ends
chrom.breaks = combined.loss.data %>% arrange(chrom, end) %>% group_by(chrom) %>% summarise(max(end))
chrom.breaks$chrom = as.integer(as.character(chrom.breaks$chrom))
chrom.breaks = chrom.breaks[order(chrom.breaks$chrom), ]
colnames(chrom.breaks)[2] = 'end'

chrom.out = combined.loss.data$bin[which(combined.loss.data$end %in% chrom.breaks$end)]


loss.plot = ggplot(combined.loss.data, aes(x = bin)) +
  geom_bar(aes(y = ratio.tumor), stat = 'identity', alpha = 0.5, width = 1, fill = tumor.color) +
  geom_bar(aes(y = ratio.csf), stat = 'identity', alpha = 1, width = 1, fill = csf.color) +
  geom_vline(xintercept = chrom.out, linetype = "dashed", colour = "grey") +
  geom_path(aes(y = ratio.tumor), group = 1, col = 'black') +
  theme_cowplot() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        panel.background = element_blank(),
        legend.position = "right",
        #aspect.ratio = 0.25,
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  scale_x_discrete(drop = FALSE) +
  scale_y_reverse(limits=c(1,0)) +
  ggtitle("") +
  ylab("Prop. of Patients\nWith Loss/LOH")


## combine everything together
genomewide = plot_grid(gain.plot, loss.plot, ncol = 1, align = 'v', axis = 'tblr', rel_heights = c(4, 4))
ggsave(filename = '05_Plots/GenomeWidePlot.pdf', plot = genomewide, device = 'pdf', width = 14, height = 9.5, bg = 'white')

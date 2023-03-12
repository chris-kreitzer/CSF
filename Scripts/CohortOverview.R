##----------------+
## Cohort overview 
## CSF/solid tumor
##----------------+


## start: 11/14/2022
## chris-kreitzer


clean()
gc()
.rs.restartR()
setup(working.path = '~/Documents/MSKCC/Subhi/CSF/')


sample_match = read.csv('Data/FINAL_samples/sample_match.txt', sep = '\t')
annotation = readxl::read_excel('Data/FINAL_samples/CSF_Lastest_07102022.xlsx')

##----------------+
## only work with QC 
## TRUE samples
##----------------+
samplesQC = sample_match[which(sample_match$fit == 'pass'), ]
samplesQC = samplesQC[!grepl(pattern = '.N0.*', x = samplesQC$sample), ]
samplesQC$path = NULL
samplesQC$remark = NULL

samples_keep = merge(samplesQC, annotation[,c('Patient ID', 'Sample ID', 'Name', 'TYPE', 'ORDER', 'Sex')],
                     by.x = 'sample', by.y = 'Sample ID', all.x = T)


##----------------+
## cohort overview; 
## PRE-work; incomplete data
##----------------+
samples_keep = samples_keep[!grepl(pattern = '\\*', x = samples_keep$ORDER), ]
re = samples_keep$sample[which(samples_keep$ORDER == 1 & samples_keep$TYPE == 'CSF')]
samples_keep = samples_keep[!samples_keep$sample %in% re, ]
re_nn = c()
for(i in unique(samples_keep$PATIENT_ID)){
  if(length(samples_keep$sample[which(samples_keep$PATIENT_ID == i)]) == 1){
    re_n = i
    re_nn = c(re_nn, re_n)
  } else next
}
samples_keep = samples_keep[!samples_keep$PATIENT_ID %in% re_nn, ]

#' order by sample amount
ordering = data.frame()
for(i in unique(samples_keep$PATIENT_ID)){
  out = data.frame(id = i,
                   n = length(samples_keep$sample[which(samples_keep$PATIENT_ID == i)]))
  ordering = rbind(ordering, out)
}
ordering = ordering[order(ordering$n, decreasing = T), ]
ordering$id = factor(ordering$id, levels = ordering$id)
samples_keep$PATIENT_ID = factor(samples_keep$PATIENT_ID, levels = ordering$id)
samples_keep$ORDER = factor(samples_keep$ORDER, levels = c(1,2,3,4,5,6,7,8))

rr = ggplot(samples_keep, aes(x = ORDER, y = PATIENT_ID, fill = TYPE)) +
  geom_point(shape = 21, size = 6, stroke = 0.3) +
  scale_x_discrete(position = 'top') +
  scale_fill_manual(values = c('TUMOR' = 'grey15',
                               'CSF' = 'white')) +
  scale_color_manual(values = c('TUMOR' = 'grey35',
                                'CSF' = 'grey35')) +
  theme_minimal() +
  theme(legend.position = 'bottom',
        axis.text = element_text(colour = 'black')) +
  labs(y = '')



##----------------+
## Fraction CNAs
##----------------+
alterations = read.csv('Data/Final/numericAlterations.txt', sep = '\t')
alterations = alterations[which(alterations$id %in% samples_keep$sample), ]
fraction_out = data.frame()
for(i in unique(alterations$id)){
  id = samples_keep$PATIENT_ID[which(samples_keep$sample == i)]
  sample = i
  n = length(alterations$gene[which(alterations$id == i & alterations$n_call %in% c(-2, -1, 1, 2))])
  out = data.frame(id = id,
                   sample = sample,
                   n = n)
  fraction_out = rbind(fraction_out, out)
}

fraction_out = rbind(fraction_out, data.frame(id = 'C-NC3N8D',
                                              sample = c('s_C_NC3N8D_T903_dZ_IM7', 's_C_NC3N8D_T904_dZ_IM7'),
                                              n = c(0, 0)))
fraction_out = merge(fraction_out, samples_keep[,c('sample', 'ORDER')],
                     by.x = 'sample', by.y = 'sample', all.x = T)
fraction_out$id = factor(fraction_out$id, levels = ordering$id)

uu = ggplot(fraction_out, aes(x = id, y = n, fill = ORDER)) +
  geom_bar(stat = 'identity', position = 'fill', color = 'black', size = 0.08) +
  scale_fill_manual(values = c('1' = 'grey75',
                               '2' = '#8c131c', 
                               '3' = '#8d322b',
                               '4' = '#dc6a4d',
                               '5' = '#e9ab91',
                               '6' = '#f7e5d9',
                               '7' = '#fff9f5',
                               '8' = 'white')) +
  coord_flip() +
  theme_minimal() +
  scale_y_continuous(position = 'right', labels = c(0, 25, 50, 75, 100), 
                     expand = c(0, 0)) +
  theme(legend.position = 'bottom',
        axis.text.y = element_blank(),
        axis.text.x = element_text(color = 'black')) +
  labs(y = '', x = '')

uu

#' plot in a grid
plot_grid(rr, uu, ncol = 2, align = 'hv', rel_widths = c(1, 0.2))


clean()
gc()

setwd('~/Documents/MSKCC/11_CSF/')
database = readxl::read_excel('00_Data/Database_Final+Feb24_ck_work_2023.xlsx', sheet = 1)

pairing = database[,c('Sample ID', 'Patient ID', 'TYPE', 'ORDER', 'Fit')]
pairing = pairing[!pairing$Fit %in% c('N/Possible', 'N/Avail'), ]
pairing = pairing[!pairing$`Patient ID` %in% names(which(table(pairing$`Patient ID`) == 1)), ]
write.table(x = pairing, file = '00_Data/SamplePairing.txt', sep = '\t', row.names = F, quote = F)

pair = readxl::read_excel('00_Data/pairing.xlsx')
pair = as.data.frame(pair)
colnames(pair)[6] = 'plot'
colnames(pair)[1] = 'PatientID'


plot = ggplot(data = as.data.frame(pair), 
              aes(x = ORDER, y = PatientID, 
                  shape = TYPE)) +
  geom_point(aes(color = plot), size = 2) +
  scale_shape_manual(values = c('TUMOR' = 16,
                                'CSF' = 17)) +
  scale_color_manual(values = c('notuse' = 'grey85',
                               'use' = 'black')) +
  scale_x_discrete(position = 'top') +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 6)) +
  labs(x = 'Order', y = '')
  

plot
ggsave(filename = 'SampleSelection.pdf', plot = plot, device = 'pdf', width = 18, height = 13)



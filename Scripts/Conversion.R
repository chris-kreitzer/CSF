setwd('~/Documents/MSKCC/11_CSF/')
xx = list.files(path = '~/Documents/MSKCC/11_CSF/07_CSF_refit/', pattern = '*.summary.txt$', recursive = T, full.names = T)

all_out = data.frame()
for(i in 1:length(xx)){
  datain = read.csv(file = xx[i], sep = '\t')
  all_out = rbind(all_out, datain)
}


library(xlsx)
dim(all_out)
write.xlsx(x = all_out, file = '~/Desktop/CSF.summary.xlsx')



csf = all_out[,c(1,11:21)]
head(csf)

csf$CDKN2A[which(csf$CDKN2A == 'Diploid')] = 0
csf$CDKN2A[which(csf$CDKN2A == 'Deletion')] = -2
csf$CDKN2A[which(csf$CDKN2A == 'Hetloss')] = -1
csf$CDKN2A[which(csf$CDKN2A == 'Gain')] = 1
csf$CDKN2A[is.na(csf$CDKN2A)] = NA
unique(csf$CDKN2A)

csf$CDK4[which(csf$CDK4 == 'Diploid')] = 0
csf$CDK4[which(csf$CDK4 == 'Deletion')] = -2
csf$CDK4[which(csf$CDK4 == 'Hetloss')] = -1
csf$CDK4[which(csf$CDK4 == 'Gain')] = 1
csf$CDK4[is.na(csf$CDK4)] = NA
csf$CDK4[which(csf$CDK4 == 'Amplification')] = 2
csf$CDK4[which(csf$CDK4 == 'Ampflification')] = 2
csf$CDK4[which(csf$CDK4 == 'CNLOH')] = -1
unique(csf$CDK4)

csf$EGFR[which(csf$EGFR == 'Diploid')] = 0
csf$EGFR[which(csf$EGFR == 'Deletion')] = -2
csf$EGFR[which(csf$EGFR == 'Hetloss')] = -1
csf$EGFR[which(csf$EGFR == 'Gain')] = 1
csf$EGFR[is.na(csf$EGFR)] = NA
csf$EGFR[which(csf$EGFR == 'Amplification')] = 2
csf$EGFR[which(csf$EGFR == 'Ampflification')] = 2
csf$EGFR[which(csf$EGFR == 'CNLOH')] = -1
unique(csf$EGFR)

csf$CDK6[which(csf$CDK6 == 'Diploid')] = 0
csf$CDK6[which(csf$CDK6 == 'Deletion')] = -2
csf$CDK6[which(csf$CDK6 == 'Hetloss')] = -1
csf$CDK6[which(csf$CDK6 == 'Gain')] = 1
csf$CDK6[is.na(csf$CDK6)] = NA
csf$CDK6[which(csf$CDK6 == 'Amplification')] = 2
csf$CDK6[which(csf$CDK6 == 'Ampflification')] = 2
csf$CDK6[which(csf$CDK6 == 'CNLOH')] = -1
unique(csf$CDK6)


csf$PTEN[which(csf$PTEN == 'Diploid')] = 0
csf$PTEN[which(csf$PTEN == 'Deletion')] = -2
csf$PTEN[which(csf$PTEN == 'Hetloss')] = -1
csf$PTEN[which(csf$PTEN == 'Gain')] = 1
csf$PTEN[is.na(csf$PTEN)] = NA
csf$PTEN[which(csf$PTEN == 'Amplification')] = 2
csf$PTEN[which(csf$PTEN == 'Ampflification')] = 2
csf$PTEN[which(csf$PTEN == 'CNLOH')] = -1
unique(csf$PTEN)


csf$KIT[which(csf$KIT == 'Diploid')] = 0
csf$KIT[which(csf$KIT == 'Deletion')] = -2
csf$KIT[which(csf$KIT == 'Hetloss')] = -1
csf$KIT[which(csf$KIT == 'Gain')] = 1
csf$KIT[is.na(csf$KIT)] = NA
csf$KIT[which(csf$KIT == 'Amplification')] = 2
csf$KIT[which(csf$KIT == 'Ampflification')] = 2
csf$KIT[which(csf$KIT == 'CNLOH')] = -1
unique(csf$KIT)


csf$MDM2[which(csf$MDM2 == 'Diploid')] = 0
csf$MDM2[which(csf$MDM2 == 'Deletion')] = -2
csf$MDM2[which(csf$MDM2 == 'Hetloss')] = -1
csf$MDM2[which(csf$MDM2 == 'Gain')] = 1
csf$MDM2[is.na(csf$MDM2)] = NA
csf$MDM2[which(csf$MDM2 == 'Amplification')] = 2
csf$MDM2[which(csf$MDM2 == 'Ampflification')] = 2
csf$MDM2[which(csf$MDM2 == 'CNLOH')] = -1
unique(csf$MDM2)


csf$MET[which(csf$MET == 'Diploid')] = 0
csf$MET[which(csf$MET == 'Deletion')] = -2
csf$MET[which(csf$MET == 'Hetloss')] = -1
csf$MET[which(csf$MET == 'Gain')] = 1
csf$MET[is.na(csf$MET)] = NA
csf$MET[which(csf$MET == 'Amplification')] = 2
csf$MET[which(csf$MET == 'Ampflification')] = 2
csf$MET[which(csf$MET == 'CNLOH')] = -1
unique(csf$MET)


csf$RB1[which(csf$RB1 == 'Diploid')] = 0
csf$RB1[which(csf$RB1 == 'Deletion')] = -2
csf$RB1[which(csf$RB1 == 'Hetloss')] = -1
csf$RB1[which(csf$RB1 == 'Gain')] = 1
csf$RB1[is.na(csf$RB1)] = NA
csf$RB1[which(csf$RB1 == 'Amplification')] = 2
csf$RB1[which(csf$RB1 == 'Ampflification')] = 2
csf$RB1[which(csf$RB1 == 'CNLOH')] = -1
unique(csf$RB1)


csf$PDGFRA[which(csf$PDGFRA == 'Diploid')] = 0
csf$PDGFRA[which(csf$PDGFRA == 'Deletion')] = -2
csf$PDGFRA[which(csf$PDGFRA == 'Hetloss')] = -1
csf$PDGFRA[which(csf$PDGFRA == 'Gain')] = 1
csf$PDGFRA[is.na(csf$PDGFRA)] = NA
csf$PDGFRA[which(csf$PDGFRA == 'Amplification')] = 2
csf$PDGFRA[which(csf$PDGFRA == 'Ampflification')] = 2
csf$PDGFRA[which(csf$PDGFRA == 'CNLOH')] = -1
unique(csf$PDGFRA)


csf$KDR[which(csf$KDR == 'Diploid')] = 0
csf$KDR[which(csf$KDR == 'Deletion')] = -2
csf$KDR[which(csf$KDR == 'Hetloss')] = -1
csf$KDR[which(csf$KDR == 'Gain')] = 1
csf$KDR[is.na(csf$KDR)] = NA
csf$KDR[which(csf$KDR == 'Amplification')] = 2
csf$KDR[which(csf$KDR == 'Ampflification')] = 2
csf$KDR[which(csf$KDR == 'CNLOH')] = -1
csf$KDR[which(csf$KDR == 'none')] = 0
unique(csf$KDR)

x = t(csf)
colnames(x) = x[1,]
x = x[-1, ]
row.names(x)
str(x)
write.table(x = x, file = '~/Documents/MSKCC/11_CSF/00_Data/CSF_binary_053023.txt', sep = '\t', quote = F, row.names = T)

library(xlsx)
write.xlsx(x = x, file = '~/Desktop/CSF_binary_summary.xlsx', row.names = T, col.names = T, showNA = F)


clean()
gc()
library(data.table)
library(colorspace)
library(copynumber)
library(CNTools)
#BiocManager::install("copynumber")
#BiocManager::install("CNTools")

#### source this and load cytoband data ####
source("~/Documents/MSKCC/MSKCC_SlackFiles/My.copy.number (1).R")

sub_seg = fread('~/Documents/MSKCC/MSKCC_SlackFiles/IGV_paired_samples (1).txt', sep = '\t',
                data.table = F, stringsAsFactors = F)
colnames(sub_seg) = c("ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean")
frequency_data_orig = sub_seg

# Read in clinical data
sub_clin_df = read.delim("~/Documents/MSKCC/MSKCC_SlackFiles/csf_clin.txt")
sub_clin_df$Sample_ID = as.character(sub_clin_df$Sample_ID)
sub_clin_df$sample_type = ifelse(sub_clin_df$CDK12 == TRUE, "DMP", "CSF")
subtype_list = c("DMP", "CSF")


# Iterate through all the subtypes and create the frequency plots
for (i in 1:length(subtype_list)){
  subtype = subtype_list[i]
  inter = sub_clin_df[sub_clin_df$sample_type == subtype, ]
  inter = unique(inter$Sample_ID)
  frequency_data = frequency_data_orig[frequency_data_orig$ID %in% inter, ]
  
  cnseg = CNSeg(frequency_data)
  rdseg = getRS(cnseg, by = "region", imput = FALSE, XY = T, what = "mean")
  frq_data = rs(rdseg)
  frq_data.annot = frq_data[,1:3]
  frq_data = frq_data[, 4:(ncol(frq_data))]
  frq_data[] <- lapply(frq_data[], function(x) as.numeric(as.character(x)))
  frq_data.annot$chrom = as.numeric(as.character(frq_data.annot$chrom))
  frq_data.annot$end = as.numeric(as.character(frq_data.annot$end))
  frq_data.annot$start = as.numeric(as.character(frq_data.annot$start))
  frq_data.annot$size = frq_data.annot$end - frq_data.annot$start
  frq_data.annot
  
  #### choose a threshold of gain/loss (usually 0.2) ####
  frq_data_calls = matrix(0, nrow = nrow(frq_data), ncol = ncol(frq_data))
  rownames(frq_data_calls) = rownames(frq_data)
  colnames(frq_data_calls) = colnames(frq_data)
  frq_data_calls[frq_data > .2] = 1
  frq_data_calls[frq_data < -.2] = -1
  frq_data_calls[1:5,1:5]
  
  #### reduce the frq_data_calls (remove consecutive segments are identical (in the same chromosome!)) ####
  chrom = unique(frq_data.annot$chrom)
  identic.idx = list()
  for(y in 1:length(chrom)) {
    frq_data_calls.chr = frq_data_calls[which(frq_data.annot$chrom == chrom[y]),]
    identic = vector(length = nrow(frq_data_calls.chr))
    for(i in  1:(nrow(frq_data_calls.chr)-1)) {
      identic[i] = !identical(frq_data_calls.chr[i,], frq_data_calls.chr[i+1,])
      identic[nrow(frq_data_calls.chr)] = !identical(frq_data_calls.chr[nrow(frq_data_calls.chr)-1,], frq_data_calls.chr[nrow(frq_data_calls.chr),])
    }
    identic.idx[[y]] = identic
  }
  identic.idx = do.call(c, identic.idx)
  frq_data_calls = frq_data_calls[identic.idx,]
  frq_data.annot = frq_data.annot[identic.idx,]
  
  #' Calculate the frequency of gains and losses
  fgain.mut = apply(frq_data_calls, MARGIN = 1, function(x) {length(which(x>0))/length(x)})
  floss.mut = apply(frq_data_calls, MARGIN = 1, function(x) {length(which(x<0))/length(x)})
  
  #' Read in cytoband
  cytoband = read.table('~/Documents/MSKCC/MSKCC_SlackFiles/cytoBand_uu.txt', header = F, sep = '\t')
  
  
  # Create pdf name
  pdf_name = paste0(subtype, "_cna_frequency_figure.pdf")
  
  pdf(file = pdf_name, width = 18, height = 6)
  print(CNA_frq_plot(frq_data_annot = frq_data.annot, 
                     fgain = fgain.mut, 
                     floss = floss.mut, 
                     col =  '#b2182b', 
                     cytoband = cytoband, 
                     ylim = c(-100, 100), 
                     title = subtype))
  dev.off()
  
}

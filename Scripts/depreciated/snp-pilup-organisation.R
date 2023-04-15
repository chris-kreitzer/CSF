##-----------------
## Organizing snp-pileup function
## for running on cluster; for CSF samples
##-----------------
##
## 09/26/2022
## chris-kreitzer


clean()
.rs.restartR()
library(stringr)
library(dplyr)

DMP_bam = read.csv('~/Documents/MSKCC/dmp-2021/Genomics/key.txt', sep = '\t', header = F)
sample_pairing = read.csv('~/Documents/MSKCC/Subhi/CSF/Data/FINAL_samples/sample_pairing.txt', sep = '\t', header = F)
research = read.csv('~/Documents/MSKCC/Subhi/CSF/Data/files.txt', sep = '\t', header = F)
research$V1 = NULL

ids = read.csv('~/Documents/MSKCC/Subhi/CSF/Data/FINAL_samples/Final_CSF_Samples.txt', sep = '\t')
#colnames(ids) = ids[1, ]
#ids = ids[-1, c('SAMPLE_ID', 'PATIENT_ID')]
#ids$PATIENT_ID = ifelse(grepl('C-', ids$PATIENT_ID), ids$PATIENT_ID, paste0('C-', ids$PATIENT_ID))


dmp_path = '/juno/res/dmpcollab/dmpshare/share/irb12_245/'
research_path = '/juno/res/ci/share/schultz/'

all_out = data.frame()
for(i in unique(ids$Patient.ID)){
  sample_subset = ids[which(ids$Patient.ID == i), ]
  for(j in 1:nrow(sample_subset)){
    id = sample_subset$Sample.ID[j]
    if(grepl('_N90', id)) next
    else {
      print(id)
      if(any(grepl(id, DMP_bam$V1))){
        normal = grep(substr(id, start = 1, stop = 9), DMP_bam$V1, value = T)
        normal = grep('-N0', normal, value = T)
        normal = ifelse(length(normal) == 2, normal[1], normal)
        normal_id = str_split_fixed(normal, pattern = ',', 2)[,1]
        normal_path = str_split_fixed(normal, pattern = ',', 3)[,2]
        normal_path = paste0(dmp_path, substr(x = normal_path, start = 0, stop = 1), '/', 
                             substr(x = normal_path, start = 2, stop = 2), '/', normal_path, '.bam')
        
        
        path = grep(id, DMP_bam$V1, value = T)
        path = str_split_fixed(path[1], pattern = ',', 3)[,2]
        path = paste0(dmp_path, substr(x = path, start = 0, stop = 1), '/', 
                      substr(x = path, start = 2, stop = 2), '/', path, '.bam')
      } else if (any(grepl(pattern = id, x = research$V2))){
        path = grep(pattern = id, research$V2, value = T)
        path = grep(pattern = '.bam$', path, value = T)
        path = paste0(research_path, substring(path, 3))
        path = ifelse(length(path) == 2, path[1], path)
        normal_id = NA
        normal_path = NA
        
      } else {
        path = NA
        normal_id = NA
        normal_path = NA
      }
      
      out = data.frame(PATIENT_ID = i,
                       sample = c(id, normal_id),
                       path = c(path, normal_path))
      all_out = rbind(all_out, out)
    }
  }
}

all_out = all_out[!with(all_out, is.na(sample) & is.na(path)), ]
all_out = all_out[!with(all_out, duplicated(sample) & duplicated(path)), ]
CSF_out = all_out %>% distinct(sample, path, .keep_all = TRUE)

##-----------------
## Research samples 
## with no matched normal
##-----------------
nu = data.frame()
for(i in unique(all_out$PATIENT_ID)){
  n = length(unique(all_out$sample[which(all_out$PATIENT_ID == i)]))
  out = data.frame(id = i,
                   n = n)
  nu = rbind(nu, out)
}

na = data.frame()
for(i in unique(ids$Patient.ID)){
  n = length(unique(ids$Sample.ID[which(ids$Patient.ID == i)]))
  out = data.frame(id = i,
                   n = n)
  na = rbind(na, out)
}

colnames(na) = c('id.o', 'n.o')
ne = merge(nu, na, by.x = 'id', by.y = 'id.o', all.x = T)
m.normal = ifelse(ne$n - ne$n.o == 0, ne$id, 'NA')
m.normal = m.normal[!m.normal == 'NA']


##-----------------
## SAMPLE pairing if 
## no DMP is available
##-----------------
m.normal
missing = data.frame()
for(i in unique(m.normal)){
  try({
    if(any(grepl(pattern = substr(x = i, start = 3, stop = nchar(i)), x = sample_pairing$V2))){
      matched = sample_pairing$V1[grepl(pattern = substr(x = i, start = 3, stop = nchar(i)), x = sample_pairing$V2)]
      matched = ifelse(length(matched) > 1, matched[1], matched)
      # matched_path = grep(pattern = matched, research$V2, value = T)
      # matched_path = ifelse(length(matched_path) > 1, matched_path[1], matched_path)
      # matched_path = grep(pattern = '.bam$', x = matched_path, value = T)
      out = data.frame(PATIENT_ID = i,
                       sample = matched,
                       path = NA)
      missing = rbind(missing, out)
      rm(matched, out)
    }
  })
}

for(i in 1:nrow(missing)){
  if(any(grepl(pattern = missing$sample[i], research$V2))){
    path = grep(pattern = missing$sample[i], research$V2, value = T)
    path = grep(pattern = '.bam$', x = path, value = T)
    path = ifelse(length(path) > 1, path[1], path)
    missing$path[i] = path
  }
}


##-----------------
CSF_match = rbind(all_out, missing)
CSF_match$path = ifelse(grepl(pattern = '^\\.', CSF_match$path),
                        paste0(research_path, substr(CSF_match$path, start = 3, stop = nchar(CSF_match$path))), 
                        CSF_match$path)


write.table(x = CSF_match, file = '~/Documents/MSKCC/Subhi/CSF/Data/FINAL_samples/sample_match.txt', sep = '\t', row.names = F, quote = F)






##-----------------
#' out


# ID_match = data.frame()
# for(i in unique(CSF_out$PATIENT_ID)){
#   if(any(grepl(pattern = 'P-0', CSF_out$sample[which(CSF_out$PATIENT_ID == i)]))) {
#     out = CSF_out[which(CSF_out$PATIENT_ID == i), ]
#   }
#   else {
#     CSF_sub = CSF_out[which(CSF_out$PATIENT_ID == i), ]
#     CSF_sample = as.character(CSF_sub[1, 'sample'])
#     normal_csf = sample_pairing$V1[grepl(pattern = CSF_sample, sample_pairing$V2)]
#     normal_csf = ifelse(length(normal_csf) > 1, normal_csf[1], normal_csf)
#     normal_csf_path = grep(pattern = normal_csf, x = research$V2, value = T)
#     normal_csf_path = grep(pattern = '.bam$', x = normal_csf_path, value = T)
#     normal_csf_path = paste0(research_path, substring(normal_csf_path, first = 3))
#     
#     out = rbind(CSF_sub, data.frame(PATIENT_ID = i,
#                                     sample = normal_csf,
#                                     path = normal_csf_path))
#   }
#   ID_match = rbind(ID_match, out)
#   rm(out)
# }
# 
# ID_match = ID_match %>% distinct(sample, .keep_all = TRUE)
# 
# write.table(x = ID_match, file = '~/Documents/MSKCC/Subhi/CSF/Data/FINAL_samples/sample_match.txt', sep = '\t', row.names = F, quote = F)

#' out


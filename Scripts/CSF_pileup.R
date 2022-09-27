snp_files = read.csv('~/CSF/sample_match.txt', sep = '\t')
snp_files = snp_files[which(snp_files$PATIENT_ID %in% c('C-50463R', 'C-VFL0PR')), ]
output_path = '~/CSF/'

snp_files_out = data.frame()
for(i in unique(snp_files$PATIENT_ID)){
  try({
    print(i)
    dir.create(path = paste0(output_path, i))
    data_sub = snp_files[which(snp_files$PATIENT_ID == i), ]
    if(nrow(data_sub) == 1) next
    else {
      for(row in 1:nrow(data_sub)){
        #print(paste0('cp ', substring(data_sub$path[row], 1, nchar(data_sub$path[row]) - 1), '* ',
        #             paste0(output_path, i, '/')))
        
        system(command = paste0('cp ', substring(data_sub$path[row], 1, nchar(data_sub$path[row]) - 1), '* ',
                                paste0(output_path, i, '/')))
      }
      
      normal = list.files(path = paste0(output_path, i), pattern = 'N0|-N.', full.names = T)
      normal = grep(pattern = '.bam$', x = normal, value = T)
      normal.bai = grep(pattern = '.bai$', x = normal, value = T)
      
      tumor = list.files(path = paste0(output_path, i), pattern = '.bam$', full.names = T)
      tumor.bai = list.files(path = paste0(output_path, i), pattern = '.bai$', full.names = T)
      tumor = tumor[!tumor %in% normal]
      
      n = length(unique(tumor))
      
      if(n == 9){
        system(command = paste0('snp-pileup -A -P50 -q15 -Q20 -r5,0 -g /juno/work/ci/resources/genomes/GRCh37/facets_snps/dbsnp_137.b37__RmDupsClean__plusPseudo50__DROP_SORT.vcf.gz ',
                                output_path, i, '/', i, '__countMatrix.dat.gz ', normal, ' ', tumor[1], ' ', tumor[2], ' ',
                                tumor[3], ' ', tumor[4], ' ', tumor[5], ' ', tumor[6], ' ', tumor[7], ' ', tumor[8], ' ', tumor[9]))
        out = data.frame(Patient_ID = i,
                         sample = c(tumor[1], tumor[2], tumor[3], tumor[4], tumor[5],
                                    tumor[6], tumor[7], tumor[8], tumor[9]),
                         pileup_file = seq(2, n + 1, 1))
        system(command = paste0('rm ', normal, ' ', tumor[1], ' ', tumor[2], ' ',
                                tumor[3], ' ', tumor[4], ' ', tumor[5], ' ', tumor[6], ' ', tumor[7], ' ', tumor[8], ' ', tumor[9]))
        system(command = paste0('rm ', normal.bai, ' ', tumor.bai[1], ' ', tumor.bai[2], ' ',
                                tumor.bai[3], ' ', tumor.bai[4], ' ', tumor.bai[5], ' ', tumor.bai[6], ' ', tumor.bai[7], ' ', tumor.bai[8], ' ', tumor.bai[9]))
        
      } else if(n == 8){
        system(command = paste0('snp-pileup -A -P50 -q15 -Q20 -r5,0 -g /juno/work/ci/resources/genomes/GRCh37/facets_snps/dbsnp_137.b37__RmDupsClean__plusPseudo50__DROP_SORT.vcf.gz ',
                                output_path, i, '/', i, '__countMatrix.dat.gz ', normal, ' ', tumor[1], ' ', tumor[2], ' ',
                                tumor[3], ' ', tumor[4], ' ', tumor[5], ' ', tumor[6], ' ', tumor[7], ' ', tumor[8]))
        out = data.frame(Patient_ID = i,
                         sample = c(tumor[1], tumor[2], tumor[3], tumor[4], tumor[5],
                                    tumor[6], tumor[7], tumor[8]),
                         pileup_file = seq(2, n+1, 1))
        system(command = paste0('rm ', normal, ' ', tumor[1], ' ', tumor[2], ' ',
                                tumor[3], ' ', tumor[4], ' ', tumor[5], ' ', tumor[6], ' ', tumor[7], ' ', tumor[8]))
        system(command = paste0('rm ', normal.bai, ' ', tumor.bai[1], ' ', tumor.bai[2], ' ',
                                tumor.bai[3], ' ', tumor.bai[4], ' ', tumor.bai[5], ' ', tumor.bai[6], ' ', tumor.bai[7], ' ', tumor.bai[8]))
        
      } else if (n == 7) {
        system(command = paste0('snp-pileup -A -P50 -q15 -Q20 -r5,0 -g /juno/work/ci/resources/genomes/GRCh37/facets_snps/dbsnp_137.b37__RmDupsClean__plusPseudo50__DROP_SORT.vcf.gz ',
                                output_path, i, '/', i, '__countMatrix.dat.gz ', normal, ' ', tumor[1], ' ', tumor[2], ' ',
                                tumor[3], ' ', tumor[4], ' ', tumor[5], ' ', tumor[6], ' ', tumor[7]))
        out = data.frame(Patient_ID = i,
                         sample = c(tumor[1], tumor[2], tumor[3], tumor[4], tumor[5],
                                    tumor[6], tumor[7]),
                         pileup_file = seq(2, n+1, 1))
        system(command = paste0('rm ', normal, ' ', tumor[1], ' ', tumor[2], ' ',
                                tumor[3], ' ', tumor[4], ' ', tumor[5], ' ', tumor[6], ' ', tumor[7]))
        system(command = paste0('rm ', normal.bai, ' ', tumor.bai[1], ' ', tumor.bai[2], ' ',
                                tumor.bai[3], ' ', tumor.bai[4], ' ', tumor.bai[5], ' ', tumor.bai[6], ' ', tumor.bai[7]))
        
      } else if (n == 6) {
        system(command = paste0('snp-pileup -A -P50 -q15 -Q20 -r5,0 -g /juno/work/ci/resources/genomes/GRCh37/facets_snps/dbsnp_137.b37__RmDupsClean__plusPseudo50__DROP_SORT.vcf.gz ',
                                output_path, i, '/', i, '__countMatrix.dat.gz ', normal, ' ', tumor[1], ' ', tumor[2], ' ',
                                tumor[3], ' ', tumor[4], ' ', tumor[5], ' ', tumor[6]))
        out = data.frame(Patient_ID = i,
                         sample = c(tumor[1], tumor[2], tumor[3], tumor[4], tumor[5],
                                    tumor[6]),
                         pileup_file = seq(2, n+1, 1))
        system(command = paste0('rm ', normal, ' ', tumor[1], ' ', tumor[2], ' ',
                                tumor[3], ' ', tumor[4], ' ', tumor[5], ' ', tumor[6]))
        system(command = paste0('rm ', normal.bai, ' ', tumor.bai[1], ' ', tumor.bai[2], ' ',
                                tumor.bai[3], ' ', tumor.bai[4], ' ', tumor.bai[5], ' ', tumor.bai[6]))
        
      } else if (n == 5) {
        system(command = paste0('snp-pileup -A -P50 -q15 -Q20 -r5,0 -g /juno/work/ci/resources/genomes/GRCh37/facets_snps/dbsnp_137.b37__RmDupsClean__plusPseudo50__DROP_SORT.vcf.gz ',
                                output_path, i, '/', i, '__countMatrix.dat.gz ', normal, ' ', tumor[1], ' ', tumor[2], ' ',
                                tumor[3], ' ', tumor[4], ' ', tumor[5]))
        out = data.frame(Patient_ID = i,
                         sample = c(tumor[1], tumor[2], tumor[3], tumor[4], tumor[5]),
                         pileup_file = seq(2, n+1, 1))
        system(command = paste0('rm ', normal, ' ', tumor[1], ' ', tumor[2], ' ',
                                tumor[3], ' ', tumor[4], ' ', tumor[5]))
        system(command = paste0('rm ', normal.bai, ' ', tumor.bai[1], ' ', tumor.bai[2], ' ',
                                tumor.bai[3], ' ', tumor.bai[4], ' ', tumor.bai[5]))
        
      } else if (n == 4) {
        system(command = paste0('snp-pileup -A -P50 -q15 -Q20 -r5,0 -g /juno/work/ci/resources/genomes/GRCh37/facets_snps/dbsnp_137.b37__RmDupsClean__plusPseudo50__DROP_SORT.vcf.gz ',
                                output_path, i, '/', i, '__countMatrix.dat.gz ', normal, ' ', tumor[1], ' ', tumor[2], ' ',
                                tumor[3], ' ', tumor[4]))
        out = data.frame(Patient_ID = i,
                         sample = c(tumor[1], tumor[2], tumor[3], tumor[4]),
                         pileup_file = seq(2, n+1, 1))
        system(command = paste0('rm ', normal, ' ', tumor[1], ' ', tumor[2], ' ',
                                tumor[3], ' ', tumor[4]))
        system(command = paste0('rm ', normal.bai, ' ', tumor.bai[1], ' ', tumor.bai[2], ' ',
                                tumor.bai[3], ' ', tumor.bai[4]))
        
      } else if (n == 3) {
        system(command = paste0('snp-pileup -A -P50 -q15 -Q20 -r5,0 -g /juno/work/ci/resources/genomes/GRCh37/facets_snps/dbsnp_137.b37__RmDupsClean__plusPseudo50__DROP_SORT.vcf.gz ',
                                output_path, i, '/', i, '__countMatrix.dat.gz ', normal, ' ', tumor[1], ' ', tumor[2], ' ',
                                tumor[3]))
        out = data.frame(Patient_ID = i,
                         sample = c(tumor[1], tumor[2], tumor[3]),
                         pileup_file = seq(2, n+1, 1))
        system(command = paste0('rm ', normal, ' ', tumor[1], ' ', tumor[2], ' ',
                                tumor[3]))
        system(command = paste0('rm ', normal.bai, ' ', tumor.bai[1], ' ', tumor.bai[2], ' ',
                                tumor.bai[3]))
        
      } else if (n == 2) {
        system(command = paste0('snp-pileup -A -P50 -q15 -Q20 -r5,0 -g /juno/work/ci/resources/genomes/GRCh37/facets_snps/dbsnp_137.b37__RmDupsClean__plusPseudo50__DROP_SORT.vcf.gz ',
                                output_path, i, '/', i, '__countMatrix.dat.gz ', normal, ' ', tumor[1], ' ', tumor[2]))
        out = data.frame(Patient_ID = i,
                         sample = c(tumor[1], tumor[2]),
                         pileup_file = seq(2, n+1, 1))
        system(command = paste0('rm ', normal, ' ', tumor[1], ' ', tumor[2]))
        system(command = paste0('rm ', normal.bai, ' ', tumor.bai[1], ' ', tumor.bai[2]))
        
      } else if (n == 1) {
        system(command = paste0('snp-pileup -A -P50 -q15 -Q20 -r5,0 -g /juno/work/ci/resources/genomes/GRCh37/facets_snps/dbsnp_137.b37__RmDupsClean__plusPseudo50__DROP_SORT.vcf.gz ',
                                output_path, i, '/', i, '__countMatrix.dat.gz ', normal, ' ', tumor[1]))
        out = data.frame(Patient_ID = i,
                         sample = c(tumor[1]),
                         pileup_file = seq(2, n+1, 1))
        system(command = paste0('rm ', normal, ' ', tumor[1]))
        system(command = paste0('rm ', normal.bai, ' ', tumor.bai[1]))
      } else next
      snp_files_out = rbind(snp_files_out, out)
    }
  })
}

write.table(snp_files_out, file = '~/CSF/snp_files_out.txt', sep = '\t', row.names = F)

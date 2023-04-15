##----------------+ 
## Classify samples as
## - FULL FACETS
## - REDUCED FACETS
## - NO FIT POSSIBLE
##----------------+


setwd('~/Documents/MSKCC/11_CSF/01_countmatrices/')

xx = list.files(path = '.', full.names = F, recursive = F)
x = readxl::read_excel('../00_Data/Database_Final+Feb23_ck_2023.xlsx', sheet = 1)

for(i in 1:nrow(x)){
  if(any(grepl(pattern = x$`Sample ID`[i], x = xx))){
    print(paste0(i, ' ok'))
  } else {
    print(x$`Sample ID`[i])
    
  }
}

for(i in 1:length(xx)){
  files = list.files(path = xx[i], pattern = '_cncf.txt$')
  if(length(files) == 0){
    print(paste0(xx[i], ' corrupt'))
  } else next
}

##-- corrupt pileups: FIT not possible
'./FROZENPOOLEDNORMAL_PITT_0214.rg.md.abra.printreads__s_C_006887_S003_d.rg.md.abra.printreads.dat.gz'
'./PN_Frozen_05500_CSF1_Zp.rg.md.abra.printreads--s_C_WXENLD_S001_d.rg.md.abra.printreads.pileup'
'./PN_FROZEN_05500_CSF2_Zp.rg.md.abra.printreads--s_C_CA6XU3_L001_d.rg.md.abra.printreads.pileup'
'./PN_FROZEN_05500_CSF2_Zp.rg.md.abra.printreads--s_C_LURHRE_S001_d.rg.md.abra.printreads.pileup'
'./PN_FROZEN_05500_CSF3_Zp.rg.md.abra.printreads--s_C_852PHC_L001_d.rg.md.abra.printreads.pileup'
'./s_C_000597_NCAS_dZ_IM5.rg.md.abra.printreads--s_C_000597_L002_d.rg.md.abra.printreads.pileup'
'./s_C_DW2V1V_N901_dZ_IM5.rg.md.abra.printreads--s_C_DW2V1V_L001_d.rg.md.abra.printreads.pileup'
'./s_C_UTH1CX_N901_dZ_IM5.rg.md.abra.printreads--s_C_UTH1CX_L001_d.rg.md.abra.printreads.pileup'




##----------------+
## Facets Status:
## - full Facets
## - reduced Facets
## - fit not possible
##----------------+
files = list.files(path = '~/Documents/MSKCC/11_CSF/01_countmatrices/', full.names = F, recursive = F)

all_df = data.frame()
for(i in 1:length(files)){
  try({
    print(files[i])
    qc = list.files(path = files[i], pattern = '_qc.txt$', full.names = T)
    fit = list.files(path = files[i], pattern = '.rds$', full.names = T)
    
    if(length(fit) != 0){
      cna_fit = readRDS(file = fit)
      qc_df = facets_fit_qc(cna_fit)
      
      ##-- full fit
      qc_full = qc_df$facets_qc
      
      ##-- rescue
      purity = qc_df$valid_purity_filter_pass
      icn_segment = qc_df$icn_allelic_state_concordance_filter_pass
      
      ##-- not possible
      homdel = qc_df$homdel_filter_pass
      waterfall = qc_df$waterfall_filter_pass
      em_cncf = qc_df$em_cncf_icn_discord_filter_pass
      contamination = qc_df$contamination_filter_pass
      subclonal = qc_df$subclonal_genome_filter_pass
      
      error = c(homdel, waterfall, em_cncf, contamination, subclonal)
      
      if(any(is.na(error))){
        error = FALSE
        qc_full = FALSE
      } else {
        error = c(homdel, waterfall, em_cncf, contamination, subclonal)
        qc_full = qc_df$facets_qc
      }
      
      
      
      if(qc_full){
        out = data.frame(id = files[i],
                         fit = 'full')
      }
      
      else if (any(!error)){
        out = data.frame(id = files[i],
                         fit = 'N/Possible')
      } 
      
      else {
        out = data.frame(id = files[i],
                         fit = 'reduced')
      }
      
      
    } else {
      out = data.frame(id = files[i],
                       fit = 'N/Avail')
    }
  })
  
  all_df = rbind(all_df, out)
}


write.table(x = all_df, file = '../00_Data/QC_summaryFirst.txt', sep = '\t', row.names = F, quote = F)
chris = readxl::read_excel('../00_Data/Database_Final+Feb23_ck_2023.xlsx', sheet = 1)
chris$id = basename(chris$FacetsCountfile)
chris$fit = NULL

xx = merge(chris, all_df, by.x = 'id', by.y = 'id', all.x = T)
xx$id = NULL
write.table(x = xx, file = '../00_Data/Database.txt', sep = '\t', row.names = F, quote = F)


#' out

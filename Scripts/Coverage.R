##----------------+
## Get raw coverage of sequencing files
##----------------+
##
## start: 03/20/2023
## chris-kreitzer



getCoverage = function(bam, 
                       target_bed, 
                       genome_assembly = "hg19"){
  
  ## read in target bed table
  target_gr = rtracklayer::import(target_bed)
  if(nchar(seqlevels(target_gr)[1]) > 3){
    seqlevels(target_gr, pruning.mode = "coarse") = c("chr1","chr2","chr3","chr4","chr5",
                                                      "chr6","chr7","chr8","chr9","chr10",
                                                      "chr11","chr12","chr13","chr14",
                                                      "chr15","chr16","chr17","chr18",
                                                      "chr19","chr20","chr21","chr22",
                                                      "chrX","chrY")
  }
  else{
    seqlevels(target_gr, pruning.mode = "coarse") = c("1","2","3","4","5",
                                                      "6","7","8","9","10",
                                                      "11","12","13","14",
                                                      "15","16","17","18",
                                                      "19","20","21","22",
                                                      "X","Y")
  }
  
  target_gr = sort(target_gr)
  # remove regions overlapped with REs
  # #target_gr = removeRE(target_gr,genome_assembly)
  # target_gr = removeGap(target_gr,genome_assembly)
  # target_gr = removeHLA(target_gr,genome_assembly)
  # target_gr = removeAQP(target_gr,genome_assembly)
  
  nRegion = length(target_gr)
  cat(paste(nRegion, "non-repeats regions from", length(seqlevels(target_gr)),
            "chromosomes in the bed file.", sep = " "))
  
  ## use helper function to calculate average coverage for each region, 
  ## in order to handle large bam files, 
  ## process 1000 regions at a time to reduce memory usage
  message("calculating depth from BAM...")
  depth = NULL
  for (i in seq(1, nRegion, 1000)){
    ## report progress
    if(i%%5000 == 1 & i > 1) cat(paste(i-1, "regions processed\n"))
    end = ifelse(i + 999 > nRegion, nRegion, i + 999)
    sub_depth = calculateSubCoverage(target_gr[i:end], bam)
    depth = c(depth, sub_depth)
  }
  ## see if number of depth equals to number of regions
  if (length(depth) == nRegion) mcols(target_gr)$depth = depth
  else stop(paste("with", nRegion, "target regions, only", 
                  length(depth), "processed. Please check your input.",
                  sep=" "))
  target_gr
}


##----------------+
## calculate average coverage 
## for each range
##----------------+
calculateSubCoverage = function(range, bam){
  ## read bam file from given ranges,
  ## filter out duplicated reads, secondary reads and unmapped reads
  ## exclude reads with mapQ==0
  param = ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE, 
                                        isSecondaryAlignment=FALSE, 
                                        isDuplicate=FALSE),
                       which=range,
                       mapqFilter=1)
  ## read alignment
  sub_alignment = GenomicAlignments::readGAlignments(bam, param = param)
  ## calculate coverage
  cov = GenomicAlignments::coverage(sub_alignment)
  cov = cov[range]
  ## return average coverage for each region
  round(mean(cov))
}


xx = getCoverage(bam = bam, target_bed = target_bed, genome_assembly = 'hg19')

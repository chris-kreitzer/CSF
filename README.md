# CSF

## Finding clinical actionable CNA from CSF samples;
- there are a lot of overlaps within the DNA alignments in CSF samples;
- meaning that two reads are supporting one Base; double count (especially for short-DNA molecules);
- snp-pileup -x --> ignore thoes reads | when supplied -x (all reads are counted)   
- -A count anomalous read pairs: https://www.biostars.org/p/91183/

# finding a way to compare DMP/CSF samples;
ratio between overlaps and non-overlaps perhaps

new test case: P-0009511-T02/T03

CSF-resources: (02/08/2023):
- sample-pileup files (countmatrices): /juno/work/ci/songt/share/nandakus/pileups
- sample-bam files (T/N): /juno/res/ci/share/schultz/GBM
- CCS: FacetsFits: /juno/work/ccs/shared/resources/impact/facets/all/

MADSEQ on CSF (ref_depth)
CDKN2A - dipLogR

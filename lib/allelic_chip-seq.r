#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library('BaalChIP')
library('dplyr')
tf = args[1]
samplesheet = args[2]
vcf = args[3]
ind = args[4]

library(rtracklayer)
gr_obj =  import("hg38-blacklist.v2.bed")
library(GenomicRanges)
blacklist_hg38 = split(gr_obj, gr_obj$name)
print(ind)
hets = c('sample'=vcf)
names(hets) = ind
print(hets)

res <- BaalChIP(samplesheet=samplesheet, hets=hets)
res <- alleleCounts(res, min_base_quality=10, min_mapq=15, verbose=FALSE)

res <- QCfilter(res, RegionsToFilter="blacklist_hg38",verbose=FALSE)

counts <- BaalChIP.get(res, 'alleleCountsPerBam')
#alleleCounts are grouped by bam_name and group_name:
names(counts)
names(counts[[ind]])

res <- mergePerGroup(res)
res <- filter1allele(res)
res <- getASB(res, Iter=5000, conf_level=0.95, cores = 4, 
              RMcorrection = TRUE, 
              RAFcorrection=TRUE)
result <- BaalChIP.report(res)
outfile = gsub("%s",tf, "%s.res")
write.csv(result, file=outfile)


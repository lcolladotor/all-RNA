library('derfinder')
library('devtools')
library('rtracklayer')

## Load counts info
counts <- read.table('/dcl01/leek/data/20sim/may202015results/cross_sample_results/counts.tsv.gz', skip = 1, stringsAsFactors = FALSE)
tmp <- read.table('/dcl01/leek/data/20sim/may202015results/cross_sample_results/counts.tsv.gz', nrows = 1, stringsAsFactors = FALSE)
colnames(counts) <- c('sample', tmp[1, 1:26], paste(as.vector(unlist(tmp[1, 27:29])), collapse = '.'), paste(as.vector(unlist(tmp[1, 30:31])), collapse = '.'))
rm(tmp)

## Process counts info
counts$totalMapped <- as.integer(sapply(strsplit(counts$total.mapped.reads, ','), '[[', 1))
save(counts, file = 'counts.Rdata')

## Chrs to work with
chrs <- paste0('chr', c(1:22, 'X', 'Y'))
#chrs <- 'chrY' ## For testing

## Load mean from bigwigs
meanRail <- lapply(chrs, function(chr) {
    x <- import(BigWigFile(file.path('/dcl01/leek/data/20sim/may202015results/coverage_bigwigs', paste0('mean.', chr, '.bw'))), as = 'RleList')
    unlist(x)
})
names(meanRail) <- chrs
save(meanRail, file = 'meanRail.Rdata')

## Locate BigWig files
bws <- rawFiles('/dcl01/leek/data/20sim/may202015results/coverage_bigwigs', samplepatt = 'samp', fileterm = NULL)
names(bws) <- gsub('.bw', '', names(bws))
bws <- bws[!grepl('unique', bws)]

## Load the coverage information without filtering
## Adjust to library size of 40 million reads
totalMap <- counts$totalMapped[ match(names(bws), counts$sample) ]
names(totalMap) <- counts$sample[ match(names(bws), counts$sample) ]
covBW <- fullCoverage(files = bws, chrs = chrs, mc.cores = 4, totalMapped = totalMap, targetSize = 40e6)

## Calculate the mean from the coverage available
system.time( meanBW <- lapply(covBW, function(x) { Reduce('+', x) / ncol(x) }) )
save(meanBW, file = 'meanBW.Rdata')
rm(covBW)

## Are there differences between the means?
diffBW <- lapply(chrs, function(chr) {
    meanBW[[chr]] - meanRail[[chr]]
})

## Overall summary of differences
lapply(diffBW, summary)
## Number of bases where methods don't agree
sapply(diffBW, function(x) { length(x[x > 0])})
## Summary of difference at bases where they don't agree
lapply(diffBW, function(x) { summary(x[ x > 0] )})



## Locate BAM files
bams <- rawFiles('/dcl01/leek/data/20sim/may202015results/alignments', samplepatt = 'bam$', fileterm = NULL)
names(bams) <- gsub('alignments.|.bam', '', names(bams))
bams <- bams[!grepl('chrM|unmapped', names(bams))]

## Load data by chr
covBAM <- lapply(chrs, function(chr) {
    bams.chr <- bams[grepl(chr, names(bams))]
    names(bams.chr) <- gsub(paste0('.', chr), '', names(bams.chr))
    ## Make sure order is correct
    bams.chr <- bams.chr[match(names(bams.chr), names(bws))]
    cov.chr <- loadCoverage(files = bams.chr, chr = chr, mc.cores = 4, totalMapped = totalMap, targetSize = 40e6)$coverage
})
names(covBAM) <- chrs

## Calculate the mean from the coverage available
system.time( meanBAM <- lapply(covBAM, function(x) { Reduce('+', x) / ncol(x) }) )
save(meanBAM, file = 'meanBAM.Rdata')
rm(covBAM)

## Is the mean from bams and from BigWigs the same?
identical(meanBAM, meanBW)


## Are there differences between the means -- vs Rail?
diffBAM <- lapply(chrs, function(chr) {
    meanBAM[[chr]] - meanRail[[chr]]
})

## Overall summary of differences vs Rail
lapply(diffBAM, summary)
## Number of bases where methods don't agree
sapply(diffBAM, function(x) { length(x[x > 0])})
## Summary of difference at bases where they don't agree
lapply(diffBAM, function(x) { summary(x[ x > 0] )})

## Done!
proc.time()
options(width = 120)
session_info()

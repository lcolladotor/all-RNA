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

## Load mean from bigwigs -- by Rail
meanRail <- lapply(chrs, function(chr) {
    x <- import(BigWigFile(file.path('/dcl01/leek/data/20sim/may202015results/coverage_bigwigs', paste0('mean.', chr, '.bw'))), as = 'RleList')
    unlist(x)
})
names(meanRail) <- chrs
save(meanRail, file = 'meanRail.Rdata')

## Load median from bigwigs -- by Rail
medianRail <- lapply(chrs, function(chr) {
    x <- import(BigWigFile(file.path('/dcl01/leek/data/20sim/may202015results/coverage_bigwigs', paste0('median.', chr, '.bw'))), as = 'RleList')
    unlist(x)
})
names(medianRail) <- chrs
save(medianRail, file = 'medianRail.Rdata')

## Choose chunks to check the median
chunksMedian <- lapply(medianRail, function(x) { 
    y <- x > 0
    ## Keep only first 100 thousand non-zeros
    cumsum(y) <= 1e5 & y
})

## Locate BigWig files
bws <- rawFiles('/dcl01/leek/data/20sim/may202015results/coverage_bigwigs', samplepatt = 'samp', fileterm = NULL)
names(bws) <- gsub('.bw', '', names(bws))
bws <- bws[!grepl('unique', bws)]

## Load the coverage information without filtering
## Adjust to library size of 40 million reads
totalMap <- counts$totalMapped[ match(names(bws), counts$sample) ]
names(totalMap) <- counts$sample[ match(names(bws), counts$sample) ]
covBW <- fullCoverage(files = bws, chrs = chrs, mc.cores = 4, totalMapped = totalMap, targetSize = 40e6)

## Calculate the mean from the coverage available -- BigWig files
system.time( meanBW <- lapply(covBW, function(x) { Reduce('+', x) / ncol(x) }) )
save(meanBW, file = 'meanBW.Rdata')

## Calculate median for some chunks
medianCalc <- function(cov.chr, med.chunk) {
    cov.mat <- as.matrix(as.data.frame(cov.chr[med.chunk, ]))
    x <- apply(cov.mat, 1, median) 
}
medianBW <- mapply(medianCalc, covBW, chunksMedian, SIMPLIFY = FALSE)
save(medianBW, file = 'medianBW.Rdata')
rm(covBW)

diffPairs <- function(rail, other, index = NULL) {
    res <- lapply(chrs, function(chr) {
        if(is.null(index)) {
            rail[[chr]] - other[[chr]]
        } else {
            rail[[chr]][index[[chr]]] - other[[chr]]
        }
        
    })
    names(res) <- chrs
    return(res)
}

summaryPairs <- function(diffInfo, tolerance = 1e-03) {
    ## Overall summary of differences
    print(lapply(diffInfo, summary))
    ## Number of bases where methods don't agree
    n.tol <- sapply(diffInfo, function(x) { length(x[abs(x) > tolerance])})
    ## Summary of difference at bases beyond the tolerance level
    if(any(n.tol)) {
        print(paste('Number of bases with difference beyond tolerance', tolerance, 'is', n.tol))
        print('Summary at bases beyond tolerance level')
        print(lapply(diffInfo, function(x) { summary(x[abs(x) > tolerance] )}))
    }
    return(invisible(NULL))
}

## Are there differences between the means?
diffBW <- diffPairs(rail = meanRail, other = meanBW)
summaryPairs(diffBW)

## Are there differences between the medians?
diffBW.med <- diffPairs(rail = medianRail, other = medianBW, chunksMedian)
summaryPairs(diffBW.med)


## Locate BAM files
bams <- rawFiles('/dcl01/leek/data/20sim/may202015results/alignments', samplepatt = 'bam$', fileterm = NULL)
names(bams) <- gsub('alignments.|.bam', '', names(bams))
bams <- bams[!grepl('chrM|unmapped', names(bams))]

## Load data by chr
covBAM <- lapply(chrs, function(chr) {
    bams.chr <- bams[grepl(paste0(chr, '$'), names(bams))]
    names(bams.chr) <- gsub(paste0('.', chr), '', names(bams.chr))
    ## Make sure order is correct
    bams.chr <- bams.chr[match(names(bams.chr), names(bws))]
    cov.chr <- loadCoverage(files = bams.chr, chr = chr, mc.cores = 4, totalMapped = totalMap, targetSize = 40e6)$coverage
})
names(covBAM) <- chrs

## Calculate the mean from the coverage available -- bam files
system.time( meanBAM <- lapply(covBAM, function(x) { Reduce('+', x) / ncol(x) }) )
save(meanBAM, file = 'meanBAM.Rdata')

## Calculate median for some chunks
medianBAM <- mapply(medianCalc, covBAM, chunksMedian, SIMPLIFY = FALSE)
save(medianBAM, file = 'medianBAM.Rdata')
rm(covBAM)

## Is the mean from bams and from BigWigs the same?
identical(meanBAM, meanBW)
if(!identical(meamBAM, meanBW)) summaryPairs(diffPairs(meanBW, meanBAM))

## What about the medians?
identical(medianBAM, medianBW)
if(!identical(medianBAM, medianBW)) summaryPairs(diffPairs(medianBW, medianBAM))


## Are there differences between the means -- vs Rail?
diffBAM <- diffPairs(rail = meanRail, other = meanBAM)
summaryPairs(diffBAM)

## Are there differences between the medians?
diffBAM.med <- diffPairs(rail = medianRail, other = medianBAM, chunksMedian)
summaryPairs(diffBAM.med)

## Done!
proc.time()
options(width = 120)
session_info()

library('GenomicRanges')
library('ggplot2')
library('scales')
library('devtools')
library('ggbio')
library('TxDb.Hsapiens.UCSC.hg19.knownGene')
library('GenomeInfoDb')
library('derfinder')

## Load data
load('/dcl01/lieber/ajaffe/derRuns/railDER/railGEU/CoverageInfo/chr22CovInfo.Rdata')
load('/dcl01/lieber/ajaffe/derRuns/railDER/railGEU/regionMatrix/regionMat-cut5-chr22.Rdata')
txdb <- keepSeqlevels(TxDb.Hsapiens.UCSC.hg19.knownGene, 'chr22')

counts <- read.table('/dcl01/lieber/ajaffe/derRuns/railDER/all_of_geuvadis_read_counts_v4.4.2015', header = TRUE, sep = '\t')

print(object.size(chr22CovInfo), units = 'Gb')
print(object.size(regionMat), units = 'Mb')


## Match the names
mapped <- counts$mapped.read.count[ 
    match(gsub('\\.', '-', colnames(chr22CovInfo$coverage)),
    counts$sample.name) ]
    
## Explore regions
regs <- regionMat$chr22$regions
regs
length(regs)
summary(width(regs))

## annotate
load("/home/epi/ajaffe/GenomicStates/GenomicState.Hsapiens.ensembl.GRCh37.p12.rda")
gs <- GenomicState.Hsapiens.ensembl.GRCh37.p12$fullGenome
ensemblAnno <- annotateRegions(regs, gs)
countTable <- ensemblAnno$countTable

### Exclusively exonic, intronic and intergenic regions
dim(countTable)
strict <- rep(NA, 3)
strict[1] <- mean(countTable[,"exon"] > 0 & countTable[,"intron"] == 0 &
    countTable[,"intergenic"] == 0) # strict exon
strict[2] <- mean(countTable[,"intron"] > 0 & countTable[,"exon"] == 0 &
    countTable[,"intergenic"] == 0) # strict intron
strict[3] <- mean(countTable[,"intergenic"] > 0 & countTable[,"exon"] == 0 &
    countTable[,"intron"] == 0) # strict intergenic
names(strict) <- c('exonic', 'intronic', 'intergenic')
strict
## Not strict
1 - sum(strict)

## Find and explore top areas
i <- order(regs$area, decreasing = TRUE)
round(regs$area[head(i, 40)], 0)
width(regs[head(i, 40)])

## Might not match longest 40 regions
table(i[seq_len(40)] %in% head(order(width(regs), decreasing = TRUE), 40))

## Longest regions (from top 40) not in top 40 areas
missing.i <- head(order(width(regs), decreasing = TRUE), 40)[ !head(order(width(regs), decreasing = TRUE), 40) %in% i[seq_len(40)] ]
width(regs[missing.i])

## Make plots for top regions
## Test:
j <- i[1]; saveLocal = TRUE

sumStat <- function(f, name, coverage) {
    stat <- eval(parse(text = paste('apply(coverage, 1,', f, ')')))
    res <- data.frame(stat = stat, name = rep(name, nrow(coverage)))
    return(res)
}

exploreTop <- function(j, saveLocal = FALSE) {
    message(paste(Sys.time(), 'plotting region', j, 'with area =', round(regs$area[j])))
    window <- resize(regs[j], round(width(regs[j]) * 1.1, 0), fix = 'center')
    range <- start(window):end(window)
    covDF <- chr22CovInfo$coverage[range, ]
    print(object.size(covDF), units = 'Mb')
    
    k <- 0
    while(k < 2) {
        p.trans <- tryCatch(autoplot(txdb, which = window, 
            names.expr = 'tx_name(gene_id)'), error = function(e) { FALSE })
        if(!is.logical(p.trans)) break
        k <- k + 1
    }
    
    
    
    for(adjust in c('Raw', 'Adjusted')) {
        message(paste(Sys.time(), 'Plot version:', adjust))
        if(adjust == 'Adjusted') {
            covDF <- DataFrame(mapply(function(x, d) x / d, covDF, mapped / 80e6))
        }
        cov <- as.matrix(as.data.frame(covDF))
        
        ## Summarize data
        df <- mapply(sumStat, c('min', 'quantile, probs = 0.05', 'quantile, probs = 0.1', 'quantile, probs = 0.25', 'mean', 'median', 'quantile, probs = 0.75', 'quantile, probs = 0.9', 'quantile, probs = 0.95', 'max'), c('min', 'p5', 'p10', 'p25', 'mean', 'median', 'p75', 'p90', 'p95', 'max'), SIMPLIFY = FALSE, MoreArgs = list(coverage = cov))
        names(df) <- NULL
        df <- lapply(df, function(x) {
            cbind(x, base = range)
        })
        df <- do.call(rbind, df)

        p <- ggplot(data = df, aes(x = base, y = stat, group = name, colour = name)) + geom_line() + ylab(paste(adjust, 'Coverage')) + scale_color_brewer(name="Statistic", palette = 'Spectral') + scale_y_continuous(trans=log2_trans())
        
        if(!is.logical(p.trans)) {
            res <- tracks(p, p.trans, heights = c(3, 1), xlim = window, title = paste('Region', j), xlab = 'Chromosome 22') + theme_bw(base_size = 18) 
        } else {
            res <- tracks(p, xlim = window, title = paste('Region', j), xlab = 'Chr22') + theme_bw(base_size = 18) 
        }
        
        if(saveLocal) pdf(file = paste0('/dcl01/lieber/ajaffe/derRuns/all-RNA/exploreGEU/exploreBar/region', j, '_', adjust, '.pdf'), width = 10)   
        #paste0('~/region', j, '_', adjust, '.pdf'), width = 10)   
        print(res)
        if(saveLocal) dev.off()
    }
    return(j)
}

pdf(file = '/dcl01/lieber/ajaffe/derRuns/all-RNA/exploreGEU/exploreBar/top40Regions_Area.pdf', width = 10)
sapply(i[seq_len(40)], exploreTop)
dev.off()

pdf(file = '/dcl01/lieber/ajaffe/derRuns/all-RNA/exploreGEU/exploreBar/top40Regions_Width_notTop40Area.pdf', width = 10)
sapply(missing.i, exploreTop)
dev.off()

## Reproducibility info
options(width = 120)
session_info()

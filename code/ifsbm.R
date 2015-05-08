# IFSBM script

###################
# Using DNAcopy
###################

# If not DNAcopy already installed, run these lines first
source("http://bioconductor.org/biocLite.R")
biocLite("DNAcopy")

# Loading the DNAcopy package in your session
require(DNAcopy)
op <- par(no.readonly = TRUE)

# Loading the supp. data:
# rename first the path according to the files location on your computer
path <- "/Users/fredcommo/Documents/myProjects/IFSBM/data"
load(file.path(path, "hg19.rda"))
load(file.path(path, "geneDB.rda"))

# Reading the data
filePath <- file.path(path, "Affy_cytoScan.cyhd.CN5.CNCHP_short.txt.gz")
foo <- readLines(filePath, n=750)
idx <- grep("ProbeSet", foo)
cnSet <- read.delim(filePath, skip=idx-1, stringsAsFactors=FALSE)
dim(cnSet)
head(cnSet)

# Filtering the SNP probes: probe Ids starting with S
cnSet <- cnSet[grep("^S", cnSet$ProbeSet),]
dim(cnSet)
head(cnSet)

# Checking the Chr names
table(cnSet$Chromosome)
cnSet$Chromosome[cnSet$Chromosome=="X"] <- 23
cnSet$Chromosome <- as.numeric(cnSet$Chromosome)
table(cnSet$Chromosome)
cnSet <- cnSet[order(cnSet$Chromosome, cnSet$Position),]

# Calculating the genomic positions
# Probes are located by their chromosome position
# Having the genomic positions will help to plot the profile.
locs <- lapply(1:24, function(chr){
    l <- cnSet$Position[cnSet$Chromosome==chr]
    l + hg19$cumlen[chr]
})
locs <- do.call(c, locs)

# A quick preview using a random subset of probes
s <- sample(1:nrow(cnSet), 1e4)
plot(locs[s], cnSet$Log2Ratio[s], cex=.2, xlab="Loc", ylab=expression(Log[2](Ratio)))
abline(h = 0, col="red")

#################################################
# Segmenting using DNAcopy and the default params
#################################################

# Constructing a DNAcopy object
LR <- cnSet$Log2Ratio
Chr <- cnSet$Chromosome
cnaObj <- CNA(LR, Chr, locs, presorted = TRUE)
cnaObj <- smooth.CNA(cnaObj)

# Segmenting
segObj <- segment(cnaObj, undo.splits = "sdundo")
sprintf("I'm an object of class %s. Great job so far!", class(segObj))

# Plotting
ptcols <- c("grey65", "grey85")
plot(segObj, xmaploc=TRUE, pt.cols=ptcols)

# Checking the Undo.SD effect
# Undo.SD is a tolerance parameter
# how many standard deviation are required to keep 2 populations distinc.

# Since the DNAcopy plot function is quite slow, let's define first our own function.
myPlot <- function(object, ptcols=c("grey65", "grey85"), n=1e4, segcols="red",...){
    
    # object has to be an object of class DNAcopy
    if(!inherits(object, "DNAcopy"))
        stop("'object' has to be of class DNAcopy.\n")
    
    dat <- object$data
    s <- sort(sample(nrow(dat), n))
    lr <- dat$Sample.1[s]
    loc <- dat$maploc[s]
    chr <- dat$chrom[s]
    st <- object$output
    cols <- lapply(unique(chr), function(ii){
        if(ii%%2==0)
            return(rep(ptcols[2], sum(chr==ii)))
        return(rep(ptcols[1], sum(chr==ii)))
    })
    cols <- do.call(c, cols)
    plot(loc, lr, ylim=range(lr, na.rm=TRUE)*1.25, col=cols,...)
    abline(h=0)
    idx <- cumsum(as.numeric(table(chr)))
    abline(v=loc[idx], lty=3)
    segments(x0=st$loc.start, x1=st$loc.end, y0=st$seg.mean, col=segcols,...)
    # You could add here some code to print the Chr nums
}

par(mfrow=c(1, 2))
for(s in c(1, 10)){
    segObj <- segment(cnaObj, undo.splits = "sdundo", undo.SD = s)
    myPlot(segObj, cex=.1, lwd=5,
           main=sprintf("I've merged pops closer than %s SDev (Undo.SD)", s))
    }
par(op) # back to default graphical params

# The alpha tolerance effect
par(mfrow=c(1, 2))
for(a in c(1e-1, 1e-10)){
    segObj <- segment(cnaObj, undo.splits = "sdundo", undo.SD = .5, alpha = a)
    myPlot(segObj, cex=.1, lwd=5,
           main=sprintf("I've kept segments distinct when\nstatistic p-value < %s (alpha)", a))
}
par(op) # back to default graphical params

# Centering the profile
require(mclust)
rLR <- runmed(LR, k=101)
rLR <- sort(rLR)
idx <- seq(1, length(rLR), len=25e3)
model <- Mclust(rLR[idx])

# A homemade function to plot the mixture model
myDplot <- function(model, nclass=200, ...){
    means <- model$parameters$mean
    props <- model$parameters$pro
    s2 <- model$parameters$variance$sigmasq
    if(length(s2)==1) s2 <- rep(s2, length(means))  

    h <- hist(model$data, nclass=nclass, plot=FALSE)
    plot(h, freq=FALSE, ylim=range(h$density)*1.25,
         xlab=expression(Log[2](Ratio)),...)
    
    md <- lapply(1:length(means), function(ii){
        m <- means[ii]; s <- sqrt(s2[ii]); p <- props[ii]
        d <- dnorm(rLR, m, s)*p
        lines(rLR, d, col=ii, lwd=4)
        text(m, max(d)+.25, labels=format(m, digits=3))
        return(max(d))
    })
    ii <- which.max(md)
    cat("Maximum peak index:", ii, ", at:", means[ii], "\n")
    
    # we return here the index of the highest peak
    return(ii)
}
maxPeak <- myDplot(model, main="The EM model")


# Final profile
choice <- model$parameters$mean[maxPeak]
LR <- LR - choice

cnaObj <- CNA(LR, Chr, locs, presorted = TRUE)
cnaObj <- smooth.CNA(cnaObj)
undo.SD <- 2
alpha <- 1e-3
segObj <- segment(cnaObj, undo.splits = "sdundo", undo.SD = undo.SD, alpha = alpha)
myPlot(segObj, cex=.1, lwd=5,
       main=sprintf("Final profile using\nUndo.SD: %s, alpha: %s", undo.SD, alpha))


#################################################
# Getting the gene list within a specifi segment
#################################################

st <- segObj$output             # to extract the segmentation table
chr <- 17                       # let's have a look on chr17
thresh <- 2                     # let's consider the segments with Log2Ratio >= 2

sgt <- st[st$chrom==chr,]
subdb <- geneDB[geneDB$chr==chr,]
ii <- which(sgt$seg.mean>=thresh)

# If length(ii)>1, this part should be put within a loop
s <- sgt$loc.start[ii]
e <- sgt$loc.end[ii]

idx <- which(s <= subdb$genomStart & subdb$genomEnd<=e)
subdb[idx, c("symbol", "fullName")]

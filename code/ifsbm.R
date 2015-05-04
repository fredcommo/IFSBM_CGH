# IFSBM script

###################
# Using DNAcopy
###################

require(DNAcopy)
op <- par(no.readonly = TRUE)

# Loading supp. data
path <- "/Users/fredcommo/Documents/myProjects/IFSBM/data"
load(file.path(path, "hg19.rda"))
load(file.path(path, "geneDB.rda"))

# Reading file
filePath <- file.path(path, "Affy_cytoScan.cyhd.CN5.CNCHP_short.txt.gz")
foo <- readLines(filePath, n=750)
idx <- grep("ProbeSet", foo)
cnSet <- read.delim(filePath, skip=idx-1, stringsAsFactors=FALSE)
dim(cnSet)
head(cnSet)

# Filtering SNP probes
cnSet <- cnSet[grep("^S", cnSet$ProbeSet),]
dim(cnSet)
head(cnSet)

# # Reducing the data
# cnSet <- cnSet[seq(1, nrow(cnSet), by=4),]
# dim(cnSet)
# head(cnSet)

# Checking Chr names
table(cnSet$Chromosome)
cnSet$Chromosome[cnSet$Chromosome=="X"] <- 23
#cnSet$Chromosome[cnSet$Chromosome=="Y"] <- 24
cnSet$Chromosome <- as.numeric(cnSet$Chromosome)
table(cnSet$Chromosome)
cnSet <- cnSet[order(cnSet$Chromosome, cnSet$Position),]

# Calculating genomic positions
locs <- lapply(1:24, function(chr){
    l <- cnSet$Position[cnSet$Chromosome==chr]
    l + hg19$cumlen[chr]
})
locs <- do.call(c, locs)
s <- sample(1:nrow(cnSet), 1e4)

# A quick preview
plot(locs[s], cnSet$Log2Ratio[s], cex=.2, xlab="Loc", ylab=expression(Log[2](Ratio)))
abline(h = 0, col="red")

# Segmenting using DNAcopy and the default params
    # Constructing a DNAcopy object
LR <- cnSet$Log2Ratio
Chr <- cnSet$Chromosome
cnaObj <- CNA(LR, Chr, locs, presorted = TRUE)
cnaObj <- smooth.CNA(cnaObj)
segObj <- segment(cnaObj, undo.splits = "sdundo")

    # Plotting
ptcols <- c("grey65", "grey85")
plot(segObj, xmaploc=TRUE, pt.cols=ptcols)

# Checking the Undo.SD effect
par(mfrow=c(1, 2))
for(s in c(1, 10)){
    segObj <- segment(cnaObj, undo.splits = "sdundo", undo.SD = s)
    plot(segObj, xmaploc=TRUE, pt.cols=ptcols)
    legend("bottomright", legend=sprintf("Undo.SD: %s", s))
}
par(op)

# The alpha tolerance effect
par(mfrow=c(1, 2))
for(a in c(1e-1, 1e-10)){
    segObj <- segment(cnaObj, undo.splits = "sdundo", undo.SD = .5, alpha = a)
    plot(segObj, xmaploc=TRUE, pt.cols=c("grey75", "grey95"))
    legend("bottomright", legend=sprintf("alpha: %s", a))
}
par(op)

# Centering the profile
require(mclust)
rLR <- runmed(LR, k=101)
rLR <- sort(rLR)
idx <- seq(1, length(rLR), len=25e3)
model <- Mclust(rLR[idx])
means <- model$parameters$mean
props <- model$parameters$pro
s2 <- model$parameters$variance$sigmasq
if(length(s2)==1) s2 <- rep(s2, length(means))  

plot(model, what="density",
     xlab=expression(Log[2](Ratio)),
     main="EM modeling")
for(ii in 1:length(means)){
    m <- means[ii]; s <- sqrt(s2[ii]); p <- props[ii]
    d <- dnorm(rLR, m, s)*p
    lines(rLR, d, col=ii, lwd=4)
    text(m, max(d)+.25, labels=format(m, digits=3))
}

# Final profile
choice <- -0.134
LR <- LR - choice

cnaObj <- CNA(LR, Chr, locs, presorted = TRUE)
cnaObj <- smooth.CNA(cnaObj)
segObj <- segment(cnaObj, undo.splits = "sdundo", undo.SD = 1, alpha = 1e-2)
plot(segObj, xmaploc=TRUE, pt.cols=ptcols)

# Getting the gene list within a specifi segment
st <- segObj$output
chr <- 17
sgt <- st[st$chrom==chr,]
subdb <- geneDB[geneDB$chr==chr,]
s <- sgt$loc.start[6]
e <- sgt$loc.end[6]

idx <- which(s <= subdb$genomStart & subdb$genomEnd<=e)
subdb[idx, c("symbol", "fullName")]

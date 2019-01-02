pkgName <- c("Daim", "BSgenome.Mmusculus.UCSC.mm10", "org.Mm.eg.db", "TxDb.Mmusculus.UCSC.mm10.knownGene")
pkgLoad <- lapply(pkgName, library, character.only = TRUE)

rawFile <- system.file("vignette", "rawData.rds", package = "Daim")
rawData <- readRDS(rawFile)

groupFactor <- factor(c(0, 0, 1, 1), levels = c(0, 1))

png("vignette/plotCoverage.png", width = 360, height = 360)
plotCoverage(rawData, group = groupFactor)
dev.off()

png("vignette/plotComplexity.png", width = 360, height = 360)
plotComplexity(rawData, group = groupFactor)
dev.off()

png("vignette/plotEnrichment.png", width = 360, height = 360)
plotEnrichment(rawData, group = groupFactor)
dev.off()

normData <- normalizeBias(rawData, group = groupFactor)

png("vignette/plotCorrelation.png", width = 480, height = 480)
par(mfrow = c(2, 2))
plotCorrelation(normData, contrast = c(1, 2))
plotCorrelation(normData, contrast = c(3, 4))
plotCorrelation(normData, contrast = c(1, 3))
plotCorrelation(normData, contrast = c(2, 4))
dev.off()

png("vignette/plotMA.png", width = 480, height = 480)
par(mfrow = c(2, 2))
plotMA(normData, contrast = c(1, 2))
plotMA(normData, contrast = c(3, 4))
plotMA(normData, contrast = c(3, 1))
plotMA(normData, contrast = c(4, 2))
dev.off()

png("vignette/plotDensity.png", width = 360, height = 360)
plotDensity(normData, group = groupFactor)
dev.off()

png("vignette/plotPCA.png", width = 360, height = 360)
plotPCA(normData, group = groupFactor)
dev.off()

png("vignette/plotMDS.png", width = 360, height = 360)
plotMDS(normData, group = groupFactor)
dev.off()

# Write Dam abundances
#writeAssay(normData, file = "vignette/Dam1.bigWig", sample = 1)
#writeAssay(normData, file = "vignette/Dam2.bigWig", sample = 2)

# Write Dam-fusion abundances
#writeAssay(normData, file = "vignette/Fusion1.bigWig", sample = 3)
#writeAssay(normData, file = "vignette/Fusion2.bigWig", sample = 4)

# Write Dam-fusion / Dam abundances
#writeRatio(normData, file = "vignette/Ratio1.bigWig", contrast = c(3, 1))
#writeRatio(normData, file = "vignette/Ratio2.bigWig", contrast = c(4, 2))

# Call DNA binding sites
bindSite <- callPeaks(normData, alpha = 0.05, lfc = log2(1.1))
#writeBroad(bindSite, "vignette/bindSite.broadPeak")

# Compute background methylation
inputData <- computeInput(normData)

# Export Input abundances
#writeAssay(inputData, "vignette/Input1.bigWig", sample = 1)
#writeAssay(inputData, "vignette/Input2.bigWig", sample = 2)

# Call DNA accessibility sites
openSite <- callPeaks(inputData, alpha = 0.05, lfc = log2(1.2))
#writeBroad(openSite, "vignette/openSite.broadPeak")

# Annotate peaks
bindSite <- annotatePeaks(bindSite, genome = "mm10")

png("vignette/plotTSS-abs.png", width = 360, height = 360)
plotTSS(bindSite, mode = "absolute")
dev.off()

png("vignette/plotTSS-rel.png", width = 360, height = 360)
plotTSS(bindSite, mode = "relative")
dev.off()

png("vignette/plotFeature.png", width = 360, height = 360)
plotFeature(bindSite)
dev.off()

# Import PWM data from JASPAR
motifPWM <- as.matrix(read.table("../backup/MA0142.1.pfm", skip = 1))
rownames(motifPWM) <- c("A", "C", "G", "T")

# Search for the MA0142.1 motif
bindSite <- searchMotif(BSgenome.Mmusculus.UCSC.mm10, ranges = bindSite, matrix = motifPWM)

# Center motif
bindSite <- centerMotif(bindSite, size = 500)

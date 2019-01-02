#' Annotate Regions in the Genome
#'
#' Annotate regions in the genome with nearest TSS and genomic feature.
#'
#' @param ranges A \code{GRanges} object of peak calls.
#' @param genome A character string specifying the genome assembly.
#' @return A \code{GRanges} object with additional metadata columns: \code{entrezTSS}, \code{symbolTSS}, \code{distanceTSS}, and \code{genomicFeature}.

annotatePeaks <- function(ranges, genome) {

    # Import genome information
    genomeFile <- system.file("data", "BiocManager.csv", package = "Daim")
    genomeInfo <- utils::read.csv(genomeFile, stringsAsFactors = FALSE)

    # Subset genome information
    genomeName <- match.arg(genome, choices = genomeInfo$assemblyName)
    genomeInfo <- genomeInfo[genomeInfo$assemblyName == genomeName, ]

    # Attach required packages
    pkgNames <- c(genomeInfo$TxDb, genomeInfo$OrgDb)
    pkgApply <- lapply(pkgNames, library, character.only = TRUE)

    # Define packages namespace
    annPkg <- eval(parse(text = genomeInfo$TxDb))
    orgPkg <- eval(parse(text = genomeInfo$OrgDb))

    # Get standard chromosome names
    matchNames <- c("M", "_")
    matchNames <- paste(matchNames, collapse = "|")
    chromNames <- seqlevels(annPkg)
    chromNames <- grep(matchNames, chromNames, invert = TRUE, value = TRUE)

    # Keep standard chromosome names
    annPkg <- GenomeInfoDb::keepSeqlevels(annPkg, chromNames)

    # Annotate nearest TSS
    ranges <- annotatePeaks.nearestTSS(ranges, genomeInfo)

    # Annotate genomic feature
    ranges <- annotatePeaks.overlapRef(ranges, genomeInfo)

    # Return annotated ranges
    ranges

}

annotatePeaks.nearestTSS <- function(ranges, genome) {

    # Define packages namespace
    annPkg <- eval(parse(text = genome$TxDb))
    orgPkg <- eval(parse(text = genome$OrgDb))

    # Define packages accessors
    columnValue <- genome$columnValue
    keyType <- genome$keyType

    # Extract TSS regions
    regionsTSS <- promoters(
        x = annPkg,
        downstream = 1000,
        upstream = 100,
        columns = columnValue,
        use.names = FALSE
    )

    # Extract ENTREZ identifiers
    entrezIds <- unlist(mcols(regionsTSS)[, columnValue])

    # Extract TSS symbols
    symbolsTSS <- AnnotationDbi::mapIds(
        x = orgPkg,
        keys = entrezIds,
        column = "SYMBOL",
        keytype = keyType
    )

    # Set TSS names to symbols
    names(regionsTSS) <- symbolsTSS[entrezIds]

    # Extract SYMBOL identifiers
    symbolIds <- names(regionsTSS)

    # Find nearest TSS
    posHits <- follow(ranges, regionsTSS)
    negHits <- precede(ranges, regionsTSS)

    # Find nearest TSS
    idxNear <- distanceToNearest(ranges, regionsTSS)
    idxHits <- subjectHits(idxNear)

    # Compute relative distance
    idxDist <- mcols(idxNear)$distance
    idxDist[posHits == idxHits] <- idxDist[posHits == idxHits]
    idxDist[negHits == idxHits] <- -idxDist[negHits == idxHits]

    # Extract hits index
    queryHits <- queryHits(idxNear)
    subjectHits <- subjectHits(idxNear)

    # Allocate results columns
    ranges$entrezTSS <- NA
    ranges$symbolTSS <- NA
    ranges$distanceTSS <- NA

    # Furnish results columns
    ranges$entrezTSS[queryHits] <- entrezIds[subjectHits]
    ranges$symbolTSS[queryHits] <- symbolIds[subjectHits]
    ranges$distanceTSS[queryHits] <- idxDist

    # Return annotated ranges
    ranges

}

annotatePeaks.overlapRef <- function(ranges, genome) {

    # Define package namespace
    annPkg <- eval(parse(text = genome$TxDb))
    orgPkg <- eval(parse(text = genome$OrgDb))

    # Extract required regions
    regionsCHR <- as(seqinfo(annPkg), "GRanges")
    regionsTXS <- transcripts(annPkg)

    # Extract annotated regions
    regionsALL <- list(
        TSS = promoters(regionsTXS, downstream = 1000, upstream = 100),
        TES = promoters(invertStrand(regionsTXS), downstream = 1000, upstream = 100),
        Exon = exons(annPkg),
        Intron = unlist(intronsByTranscript(annPkg)),
        UTR5 = unlist(fiveUTRsByTranscript(annPkg)),
        UTR3 = unlist(threeUTRsByTranscript(annPkg)),
        Intergenic = GenomicRanges::setdiff(regionsCHR, reduce(regionsTXS), ignore.strand = TRUE)
    )

    # Clean annotated regions
    cleanRange <- function(x) unname(granges(x))
    regionsALL <- lapply(regionsALL, cleanRange)
    regionsALL <- sort(unlist(GRangesList(regionsALL)))

    # Find overlapping annotations
    overlapALL <- findOverlaps(ranges, regionsALL, ignore.strand = TRUE, select = "all")
    queryHits <- queryHits(overlapALL)
    subjectHits <- subjectHits(overlapALL)

    # Compute overlapping fraction
    hitsInter <- pintersect(ranges[queryHits], regionsALL[subjectHits])
    fracInter <- width(hitsInter) / width(ranges[queryHits])
    fracInter <- split(fracInter, queryHits)
    keepInter <- unlist(lapply(fracInter, function(x) x == max(x)))

    # Filter for maximum overlap
    overlapALL <- overlapALL[keepInter]
    queryHits <- queryHits(overlapALL)
    subjectHits <- subjectHits(overlapALL)

    # Select overlap by feature
    hitsOrder <- c("TSS", "TES", "Exon", "UTR5", "UTR3", "Intron", "Intergenic")
    subjectNames <- names(regionsALL)[subjectHits]
    subjectNames <- split(subjectNames, queryHits)
    subjectNames <- sapply(subjectNames, function(x) head(x[order(match(x, hitsOrder))], 1))
    subjectIndex <- as.numeric(names(subjectNames))

    # Allocate results column
    ranges$genomicFeature <- NA

    # Furnish results column
    ranges$genomicFeature[subjectIndex] <- subjectNames

    # Return annotated ranges
    ranges

}

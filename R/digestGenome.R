#' Restriction digest of genome
#'
#' Extract DpnII restriction fragments from a genome.
#'
#' @param object A \code{BSgenome}, \code{DNAStringSet}, or \code{FaFile} object.
#'
#' @return A \code{\link[GenomicRanges]{GRanges}} object of DpnII restriction fragments.

digestGenome <- function(object) {

    # Check argument missing
    if (missing(object)) {
        stop("`object` is missing, with no default.", call. = FALSE)
    }

    # Check argument class
    validClass <- c("BSgenome", "DNAStringSet", "FaFile")
    wrongClass <- !(class(object) %in% validClass)
    if (wrongClass) {
        text <- paste("*", validClass, collapse = "\n")
        stop("`object` must be of class:\n", text, call. = FALSE)
    }

    # Dispatch relevant method
    UseMethod("digestGenome", object)
}

digestGenome.BSgenome <- function(object) {

    # Get genome information
    organismName <- organism(object)
    providerName <- provider(object)

    # Get standard chromosome names
    matchNames <- c("M", "_")
    matchNames <- paste(matchNames, collapse = "|")
    chromNames <- seqlevels(object)
    chromNames <- grep(matchNames, chromNames, invert = TRUE, value = TRUE)

    # Get major chromosome sizes
    chromSizes <- seqlengths(object)[chromNames]

    # Create list to collect ranges
    chromCount <- length(chromNames)
    rangesList <- vector("list", chromCount)
    rangesList <- setNames(rangesList, chromNames)

    # Create list to collect metadata
    mcolsList <- vector("list", chromCount)
    mcolsList <- setNames(rangesList, chromNames)

    # Append fragment ranges to list
    for (chromName in chromNames) {

        # Find motif match
        motifPattern <- DNAString("GATC")
        motifSubject <- object[[chromName]]
        motifMatches <- matchPattern(motifPattern, motifSubject)

        # Create fragment ranges
        matchStarts <- c(1L, start(motifMatches) + 2L)
        matchEnds <- c(start(motifMatches) + 1L, chromSizes[[chromName]])
        matchRanges <- GRanges(chromName, IRanges(matchStarts, matchEnds))

        # Append ranges to list
        rangesList[[chromName]] <- matchRanges

        # Profile fragment ranges
        matchViews <- BSgenomeViews(object, matchRanges)
        matchProbs <- letterFrequency(matchViews, letters = "GC", as.prob = TRUE)
        matchProbs <- matchProbs[, "G|C"]
        matchWidth <- width(matchViews)

        # Append profile to list
        mcolsList[[chromName]] <- DataFrame(digestBySize = matchWidth, digestByProb = matchProbs)

    }

    # Unlist fragment ranges
    rangesList <- GRangesList(rangesList)
    rangesData <- unlist(rangesList, use.names = FALSE)

    # Add metadata to ranges
    mcols(rangesData) <- do.call(rbind, mcolsList)

    # Add genome information to ranges
    chromInfo <- seqinfo(object)[chromNames]
    seqinfo(rangesData) <- chromInfo

    # Return fragment ranges
    rangesData
}

digestGenome.DNAStringSet <- function(object) {

    # Get standard chromosome names
    matchNames <- c("M", "_")
    matchNames <- paste(matchNames, collapse = "|")
    chromNames <- names(object)
    chromNames <- grep(matchNames, chromNames, invert = TRUE, value = TRUE)

    # Get major chromosome sizes
    chromSizes <- seqlengths(object)[chromNames]

    # Create list to collect ranges
    chromCount <- length(chromNames)
    rangesList <- vector("list", chromCount)
    rangesList <- setNames(rangesList, chromNames)

    # Create list to collect metadata
    mcolsList <- vector("list", chromCount)
    mcolsList <- setNames(rangesList, chromNames)

    # Append fragment ranges to list
    for (chromName in chromNames) {

        # Find motif match
        motifPattern <- DNAString("GATC")
        motifSubject <- object[[chromName]]
        motifMatches <- matchPattern(motifPattern, motifSubject)

        # Create fragment ranges
        matchStarts <- c(1L, start(motifMatches) + 2L)
        matchEnds <- c(start(motifMatches) + 1L, chromSizes[[chromName]])
        matchRanges <- GRanges(chromName, IRanges(matchStarts, matchEnds))

        # Append ranges to list
        rangesList[[chromName]] <- matchRanges

        # Profile fragment ranges
        matchViews <- object[matchRanges]
        matchProbs <- letterFrequency(matchViews, letters = "GC", as.prob = TRUE)
        matchProbs <- matchProbs[, "G|C"]
        matchWidth <- width(matchViews)

        # Append profile to list
        mcolsList[[chromName]] <- DataFrame(digestBySize = matchWidth, digestByProb = matchProbs)

    }

    # Unlist fragment ranges
    rangesList <- GRangesList(rangesList)
    rangesData <- unlist(rangesList, use.names = FALSE)

    # Add metadata to ranges
    mcols(rangesData) <- do.call(rbind, mcolsList)

    # Add genome information to ranges
    chromInfo <- seqinfo(object)[chromNames]
    seqinfo(rangesData) <- chromInfo

    # Return fragment ranges
    rangesData
}

digestGenome.FaFile <- function(object) {

    # Open fasta and index files
    fastaFile <- open(object)
    indexFile <- scanFaIndex(fastaFile)

    # Assign a relevant genome identifier
    fastaPath <- path(fastaFile)
    genome(indexFile) <- tools::file_path_sans_ext(fastaPath, compression = TRUE)

    # Read sequence names and widths
    chromInfo <- seqinfo(indexFile)

    # Get standard chromosome names
    matchNames <- c("M", "_")
    matchNames <- paste(matchNames, collapse = "|")
    chromNames <- seqnames(chromInfo)
    names(indexFile) <- chromNames
    chromNames <- grep(matchNames, chromNames, invert = TRUE, value = TRUE)

    # Get major chromosome sizes
    chromSizes <- seqlengths(chromInfo)[chromNames]

    # Create list to collect ranges
    chromCount <- length(chromNames)
    rangesList <- vector("list", chromCount)
    rangesList <- setNames(rangesList, chromNames)

    # Create list to collect metadata
    mcolsList <- vector("list", chromCount)
    mcolsList <- setNames(rangesList, chromNames)

    # Append fragment ranges to list
    for (chromName in chromNames) {

        # Find motif match
        motifPattern <- DNAString("GATC")
        motifSubject <- scanFa(fastaFile, param = indexFile[chromName])[[1]]
        motifMatches <- matchPattern(motifPattern, motifSubject)

        # Create fragment ranges
        matchStarts <- c(1L, start(motifMatches) + 2L)
        matchEnds <- c(start(motifMatches) + 1L, chromSizes[[chromName]])
        matchRanges <- GRanges(chromName, IRanges(matchStarts, matchEnds))

        # Append ranges to list
        rangesList[[chromName]] <- matchRanges

        # Profile fragment ranges
        matchViews <- scanFa(fastaFile, matchRanges)
        matchProbs <- letterFrequency(matchViews, letters = "GC", as.prob = TRUE)
        matchProbs <- matchProbs[, "G|C"]
        matchWidth <- width(matchViews)

        # Append profile to list
        mcolsList[[chromName]] <- DataFrame(digestBySize = matchWidth, digestByProb = matchProbs)

    }

    # Unlist fragment ranges
    rangesList <- GRangesList(rangesList)
    rangesData <- unlist(rangesList, use.names = FALSE)

    # Add metadata to ranges
    mcols(rangesData) <- do.call(rbind, mcolsList)

    # Add genome information to ranges
    chromInfo <- seqinfo(object)[chromNames]
    seqinfo(rangesData) <- chromInfo

    # Return fragment ranges
    rangesData
}

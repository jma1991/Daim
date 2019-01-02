#' Known motif search
#'
#' Search peak regions for known motifs using a position weight matrix.
#'
#' @param object A \code{BSgenome}, \code{DNAStringSet}, or \code{FaFile} object.
#' @param ranges A \code{GRanges} object of peak regions.
#' @param matrix A position weight \code{matrix} of the known motif.
#' @return A \code{GRanges} object with additional metadata columns:
#' \item{motifCenter}{The position of the motif center.}
#' \item{motifScore}{The score of the motif match.}

searchMotif <- function(object, ranges, matrix) {

    # Check argument missing
    if (missing(object)) {
        stop("`object` is missing, with no default.", call. = FALSE)
    }
    if (missing(ranges)) {
        stop("`ranges` is missing, with no default.", call. = FALSE)
    }
    if (missing(matrix)) {
        stop("`matrix` is missing, with no default.", call. = FALSE)
    }

    # Check argument class
    if (!is(ranges, "GRanges")) {
        stop("`object` must be of class `GRanges`.", call. = FALSE)
    }
    if (!is(matrix, "matrix")) {
        stop("`matrix` must be of class `matrix`.", call. = FALSE)
    }

    # Check argument value
    if (!all(rownames(matrix) == c("A", "C", "G", "T"))) {
        stop('rownames of `matrix` must equal c("A" "C" "G" "T")', call. = FALSE)
    }

    # Dispatch relevent method
    UseMethod("searchMotif", object)
}

searchMotif.BSgenome <- function(object, ranges, matrix) {

    # Get range sequences
    rangesSeq <- getSeq(object, ranges)

    # Join into a master sequence
    masterSeq <- unlist(rangesSeq)
    masterLen <- length(masterSeq)

    # Create range views on master sequence
    rangesIndex <- seq_len(length(ranges) - 1)
    rangesWidth <- cumsum(width(ranges))
    rangesStart <- c(1L, rangesWidth[rangesIndex] + 1)
    rangesEnd <- c(rangesWidth[rangesIndex], masterLen)
    masterViews <- Views(masterSeq, start = rangesStart, end = rangesEnd)

    # Match PWM on master string views
    matchObject <- suppressWarnings(matchPWM(matrix, masterViews, min.score = "50%", with.score = TRUE))

    # Return early if no matches
    if (length(matchObject) == 0) {
        message("No matches found!")
        return(invisible())
    }

    # Find motif matches on seperate ranges
    rangesMatch <- IRanges(start = start(matchObject), end = end(matchObject))
    hitsObject <- findOverlaps(rangesMatch, masterViews)

    # Find highest score for each range
    scoreData <- data.frame(query = queryHits(hitsObject), score = mcols(matchObject)$score)
    scoreList <- split(scoreData, subjectHits(hitsObject))
    keepScore <- sapply(scoreList, function(x) x$query[which.max(x$score)])

    # Select range with highest score
    matchObject <- matchObject[keepScore]
    rangesMatch <- rangesMatch[keepScore]
    hitsObject <- hitsObject[keepScore]

    # Extract motif and range index
    matchIndex <- queryHits(hitsObject)
    viewsIndex <- subjectHits(hitsObject)

    # Extract start position
    matchStart <- start(matchObject)
    rangeStart <- start(ranges)[viewsIndex]
    viewsStart <- start(masterViews)[viewsIndex]

    # Calculate motif position
    motifStart <- rangeStart + (matchStart - viewsStart)
    motifEnd <- motifStart + width(rangesMatch)
    motifCenter <- motifStart + round(width(rangesMatch) / 2)
    motifScore <- mcols(matchObject)$score

    # Append motifs to metadata
    mcolsData <- mcols(ranges)
    dummyData <- rep(NA, length(ranges))
    motifData <- DataFrame(motifCenter = dummyData, motifScore = dummyData)
    motifData$motifCenter[viewsIndex] <- motifCenter
    motifData$motifScore[viewsIndex] <- motifScore
    mcols(ranges) <- cbind(mcolsData, motifData)

    # Return annotated object
    ranges

}

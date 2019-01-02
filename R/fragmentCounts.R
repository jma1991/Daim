#' Count reads in restriction fragments
#'
#' Assign mapped sequencing reads to restriction fragments.
#'
#' @param reads A character vector of paths to BAM files.
#' @param fragments A \code{GRanges} object of restriction fragments.
#' @return A \code{RangedSummarizedExperiment} object.

fragmentCounts <- function(reads, fragments) {

    # Check argument missing
    if (missing(reads)) {
        stop("`reads` is missing, with no default.", call. = FALSE)
    }
    if (missing(fragments)) {
        stop("`fragments` is missing, with no default.", call. = FALSE)
    }

    # Check argument class
    if (!is(reads, "character")) {
        stop("`reads` must be a character vector.", call. = FALSE)
    }
    if (!is(fragments, "GRanges")) {
        stop("`fragments` must be a GRanges object.", call. = FALSE)
    }

    # Check BAM files exist
    bamExists <- file.exists(reads)
    if (!all(bamExists)) {
        text <- paste("*", reads[!bamExists], collapse = "\n")
        stop("Each file in `reads` must exist:\n", text, call. = FALSE)
    }

    # Index required BAM files
    baiPaths <- paste0(reads, ".bai")
    baiExists <- file.exists(baiPaths)
    if (!all(baiExists)) {
        indexBam(reads[!baiExists])
    }

    # Extract fragment ranges
    featureRanges <- granges(fragments)
    mcols(featureRanges)$id <- seq_along(featureRanges)

    # Count reads into fragments
    annotFile <- Rsubread::createAnnotationFile(featureRanges)
    invisible(capture.output(featureCounts <- Rsubread::featureCounts(
        files = reads,
        annot.ext = annotFile,
        read2pos = 5,
        ignoreDup = TRUE
    )))

    # Create experiment data
    assaysData <- list(countsData = unname(featureCounts$counts))
    sampleData <- DataFrame(bamPath = reads)

    # Combine experiment data
    exptData <- SummarizedExperiment(
        assays = assaysData,
        rowRanges = fragments,
        colData = sampleData
    )

    # Return experiment data
    return(exptData)
}

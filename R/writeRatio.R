#' Write ratio to a bigWig file
#'
#' Export ratio values over ranges to a bigWig file format.
#'
#' @param object A \code{RangedSummarizedExperiment} object.
#' @param file A character string naming a file open for writing.
#' @param contrast A character vector with exactly two elements: the name or index of the numerator level for the fold change, and the name or index of the denominator level for the change.

writeRatio <- function(object, file = "", contrast = NULL) {

    # Check argument missing
    if (missing(object)) {
        stop("`object` is missing, with no default.", call. = FALSE)
    }
    if (missing(file)) {
        stop("`file` is missing, with no default.", call. = FALSE)
    }
    if (missing(contrast)) {
        stop("`contrast` is missing, with no default.", call. = FALSE)
    }

    # Check argument class
    if (class(object) != "RangedSummarizedExperiment") {
        stop("`object` must be of class RangedSummarizedExperiment.", call. = FALSE)
    }
    if (class(file) != "character") {
        stop("`file` must be a single character string.", call. = FALSE)
    }

    # Check argument value
    colNames <- c(seq_len(ncol(object)), colnames(object))
    if (!(contrast[1] %in% colNames)) {
        stop("Can't find column `", contrast[1], "` in `object`.", call. = FALSE)
    }
    if (!(contrast[2] %in% colNames)) {
        stop("Can't find column `", contrast[2], "` in `object`.", call. = FALSE)
    }

    # Dispatch relevant method
    UseMethod("writeRatio", object)

}

writeRatio.RangedSummarizedExperiment <- function(object, file = "", contrast = NULL) {

    # Create fold change
    baseMean1 <- 2^assay(object, "qsmoothData")[, contrast[1]]
    baseMean0 <- 2^assay(object, "qsmoothData")[, contrast[2]]
    baseRatio <- baseMean1 - baseMean0

    # Create coverage track
    rowRanges <- granges(rowRanges(object))
    score(rowRanges) <- baseRatio

    # Export coverage track
    rtracklayer::export(rowRanges, con = file, format = "bigWig")
}

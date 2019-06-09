#' Write assay to a bigWig file
#'
#' Export assay values over ranges to a bigWig file format.
#'
#' @param object A \code{RangedSummarizedExperiment} object.
#' @param file A character string naming a file open for writing.
#' @param sample Either a column name or index.

writeAssay <- function(object, file = "", sample = NULL) {

    # Check argument missing
    if (missing(object)) {
        stop("`object` is missing, with no default.", call. = FALSE)
    }
    if (missing(file)) {
        stop("`file` is missing, with no default.", call. = FALSE)
    }
    if (missing(sample)) {
        stop("`sample` is missing, with no default.", call. = FALSE)
    }

    # Check argument class
    if (!is(object, "RangedSummarizedExperiment")) {
        stop("`object` must be of class RangedSummarizedExperiment.", call. = FALSE)
    }
    if (!is(file, "character")) {
        stop("`file` must be a single character string.", call. = FALSE)
    }
    if (!(class(sample) %in% c("character", "numeric"))) {
        stop("`sample` must be a numeric or character class.", call. = FALSE)
    }

    # Check argument value
    if (!(sample %in% c(seq_len(ncol(object)), colnames(object)))) {
        stop("Sample '", sample, "' not found.", call. = FALSE)
    }

    # Dispatch relevant method
    UseMethod("writeAssay", object)

}

writeAssay.RangedSummarizedExperiment <- function(object, file = "", sample = NULL) {

    assayNames <- assayNames(object)

    rowRanges <- rowRanges(object)

    if (assayNames[1] == "countsData") {
        score(rowRanges) <- assay(object)[, sample]
    } else if (assayNames[1] == "qsmoothData") {
        score(rowRanges) <- 2^assay(object)[, sample]
    }

    # Export coverage track
    rtracklayer::export(rowRanges, con = file, format = "bigWig")
}

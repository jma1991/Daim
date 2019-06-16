#' Write assay to a bigWig file
#'
#' Export assay values over ranges to a bigWig file format.
#'
#' @param object A \code{RangedSummarizedExperiment} object.
#' @param assay Name of assay slot.
#' @param file A character string naming a file open for writing.
#' @param sample Either a column name or index.

writeAssay <- function(object, file = "", assay = "", sample = "") {

    # Check argument missing
    if (missing(object)) {
        stop("`object` is missing, with no default.", call. = FALSE)
    }
    if (missing(file)) {
        stop("`file` is missing, with no default.", call. = FALSE)
    }
    if (missing(assay)) {
        stop("`assay` is missing, with no default.", call. = FALSE)
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
    if (!is(assay, "character")) {
        stop("`assay` must be a single character string.", call. = FALSE)
    }
    if (!(class(sample) %in% c("character", "numeric"))) {
        stop("`sample` must be a numeric or character class.", call. = FALSE)
    }

    # Check argument value
    if (!(assay %in% assayNames(object))) {
        stop("Assay '", assay, "' not found.", call. = FALSE)
    }
    if (!(sample %in% c(seq_len(ncol(object)), colnames(object)))) {
        stop("Sample '", sample, "' not found.", call. = FALSE)
    }

    # Dispatch relevant method
    if (assay == "countsData") {
        writeAssay.counts(object, file, sample)
    } else if (assay == "qsmoothData") {
        writeAssay.qsmooth(object, file, sample)
    } else {
        stop("Error!")
    }

}

writeAssay.counts <- function(object, file = "", sample = NULL) {

    # Extract assay and range data
    assayData <- assay(object, "countsData")
    rangeData <- rowRanges(object)

    # Filter assay data to chosen sample
    colIndex <- sample
    colAssay <- assayData[, colIndex]

    # Assign assay to score column
    score(rangeData) <- colAssay

    # Export coverage track
    rtracklayer::export(rowRanges, con = file, format = "bigWig")

}

writeAssay.qsmooth <- function(object, file = "", sample = NULL) {

    # Extract assay and range data
    assayData <- assay(object, "qsmoothData")
    rangeData <- rowRanges(object)

    # Filter assay data to chosen sample
    colIndex <- sample
    colAssay <- assayData[, colIndex]

    # Reverse log transformation from qsmooth
    colAssay <- 2^colAssay

    # Assign assay to score column
    score(rangeData) <- colAssay

    # Export coverage track
    rtracklayer::export(rowRanges, con = file, format = "bigWig")

}

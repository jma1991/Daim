#' Write a GRanges object as a broadPeak file
#'
#' \code{writePeak} prints a GRanges object in broadPeak format to a file.
#'
#' @param object A GRanges object to be written in broadPeak format.
#' @param file A character string naming a file open for writing.

writeBroad <- function(object, file = "") {

    # Check argument missing
    if (missing(object)) {
        stop("`object` is missing, with no default.", call. = FALSE)
    }
    if (missing(file)) {
        stop("`file` is missing, with no default.", call. = FALSE)
    }

    # Check argument class
    if (class(object) != "GRanges") {
        stop("`object` must be a GRanges object.", call. = FALSE)
    }
    if (class(file) != "character") {
        stop("`file` must be a character string.", call. = FALSE)
    }

    # Check accessors missing
    if (!("log2FoldChange" %in% names(mcols(object)))) {
        stop("Column `log2FoldChange` must exist in mcols of `object`.", call. = FALSE)
    }
    if (!("pval" %in% names(mcols(object)))) {
        stop("Column `pval` must exist in mcols of `object`.", call. = FALSE)
    }
    if (!("padj" %in% names(mcols(object)))) {
        stop("Column `padj` must exist in mcols of `object`.", call. = FALSE)
    }

    # Check accessors class
    if (!(is.numeric(mcols(object)$log2FoldChange))) {
        stop("Column `log2FoldChange` must be a numeric vector.", call. = FALSE)
    }
    if (!(is.numeric(mcols(object)$pval))) {
        stop("Column `pval` must be a numeric vector.", call. = FALSE)
    }
    if (!(is.numeric(mcols(object)$padj))) {
        stop("Column `padj` must be a numeric vector.", call. = FALSE)
    }

    # Dipsatch relevant method
    UseMethod("writeBroad", object)
}

writeBroad.GRanges <- function(object, file = "") {

    # Sanitize file path
    filePath <- file.path(file)
    filePath <- normalizePath(filePath, mustWork = FALSE)

    # Convert GRanges to broadPeak
    broadFormat <- data.frame(
        chrom = seqnames(object),
        chromStart = start(object) - 1,
        chromEnd = end(object),
        name = ".",
        score = pmin(1000, as.integer(-125 * log2(object$padj + 1e-12))),
        strand = ".",
        signalValue = object$log2FoldChange,
        pValue = -log10(object$pval),
        qValue = -log10(object$padj)
    )

    # Format scientific values
    broadFormat <- apply(broadFormat, 2, format, scientific = FALSE)
    broadFormat <- trimws(broadFormat)

    # Write to disk
    write.table(
        x = broadFormat,
        file = filePath,
        quote = FALSE,
        sep = "\t",
        row.names = FALSE,
        col.names = FALSE
    )
}

#' Center regions by motif
#'
#' Center peak regions by the position of the known motif.
#'
#' @param object A \code{GRanges} object of peak regions.
#' @return A \code{GRanges} object of peak regions centered by motif position.
#'

centerMotif <- function(object, size = NULL) {

    # Check argument missing
    if (missing(object)) {
        stop("`object` is missing, with no default", call. = FALSE)
    }

    # Check argument class
    if (!is(object, "GRanges")) {
        stop("`object` must be of class `GRanges`.", call. = FALSE)
    }

    # Check accessor missing
    if (is.null(object$motifCenter)) {
        stop("Column `motifCenter` must exist in metadata columns of `object`.", call. = FALSE)
    }

    # Dispatch relevant method
    UseMethod("centerMotif", object)

}

centerMotif.GRanges <- function(object, size = NULL) {

    # Extract regions with motif
    keepRange <- !is.na(object$motifCenter)

    # Extract region sizes
    if (is.null(size)) {
        rangeWidth <- width(object)[keepRange]
    }
    else {
        rangeWidth <- rep(size, length(keepRange))
    }

    # Replace start index with motif index
    start(object)[keepRange] <- object$motifCenter[keepRange]

    # Resize regions with motif
    object[keepRange] <- resize(object[keepRange], rangeWidth, fix = "start")

    # Return modified object
    object
}

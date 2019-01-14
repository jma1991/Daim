#' Filter Ranges By Size
#'
#' @param object A GRanges object.
#' @param min Minimum size required.
#' @param max Maximum size required.
#'

filterBySize <- function(object, min = 100, max = 100000) {

    # Check argument missing
    if (missing(object)) {
        stop("`object` is missing, with no default.", call. = FALSE)
    }

    # Check argument class
    if (!is(object, "GRanges")) {
        stop("`object` must be a GRanges object.", call. = FALSE)
    }
    if(!is.numeric(min)) {
        stop("`min` must be a numeric value", call. = FALSE)
    }
    if (!is.numeric(max)) {
        stop("`max` must be a numeric value.", call. = FALSE)
    }

    # Check argument value
    if (length(min) > 1) {
        stop("`min` must have length 1, not length ", length(min), ".", call. = FALSE)
    }
    if (length(max) > 1) {
        stop("`max` must have length 1, not length ", length(min), ".", call. = FALSE)
    }
    if (min < 0) {
        stop("`min` must be a positive number.", call. = FALSE)
    }
    if (max < 0) {
        stop("`max` must be a positive number.", call. = FALSE)
    }
    if (max < min) {
        stop("`max` must be greater than `min`.", call. = FALSE)
    }

    # Dispatch relevant method
    UseMethod("filterBySize", object)

}

filterBySize.GRanges <- function(object, min = 100, max = 100000) {

    # Extract range sizes
    rangeSize <- width(object)

    # Compute range filter
    rangeKeep <- (rangeSize > min) & (rangeSize < max)

    # Return range filter
    rangeKeep

}

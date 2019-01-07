#' Compute Methylation for Input
#'
#' Compute input restriction fragment methylation values for each Dam library.
#'
#' @param object A RangedSummarizedExperiment object.
#' @return A RangedSummarizedExperiment object.

computeInput <- function(object) {

    # Check argument missing
    if (missing(object)) {
        stop("`object` is missing, with no default.", call. = FALSE)
    }

    # Check argument class
    if (!is(object, "RangedSummarizedExperiment")) {
        stop("`object` must be a RangedSummarizedExperiment object.", call. = FALSE)
    }

    # Check argument value
    if (!"qsmoothData" %in% assayNames(object)) {
        stop("Matrix `qsmoothData` must exist in assays slot of `object`.", call. = FALSE)
    }

    # Check accessor missing
    if (is.null(object$groupFactor)) {
        stop("Column `groupFactor` must exist in colData of `object`.", call. = FALSE)
    }

    # Check accessor class
    if (!is(object$groupFactor, "factor")) {
        stop("Column `groupFactor` must be a factor.", call. = FALSE)
    }

    # Check accessor value
    if (!setequal(levels(object$groupFactor), c("0", "1"))) {
        stop("Column `groupFactor` must have levels '0' and '1' only.", call. = FALSE)
    }

    # Dispatch relevant method
    UseMethod("computeInput", object)
}

computeInput.RangedSummarizedExperiment <- function(object) {

    # Subset experiment data by control
    keepSamples <- which(object$groupFactor == 0)
    object <- object[, keepSamples]

    # Compute average count for each library
    colMeans <- colMeans(assay(object, "qsmoothData"))
    aveValue <- matrix(colMeans, nrow(object), ncol(object), byrow = TRUE)

    # Calculate background from Dam methylation
    assayData <- pmax(
        meanAve = aveValue,
        mean05k = rollAssay(object, size = 5e3),
        mean10k = rollAssay(object, size = 1e4)
    )

    # Replace assay data
    assaysData <- cbind(assayData, assay(object))
    assaysList <- list(qsmoothData = assaysData)

    # Define group columns
    groupValues <- rep(c(0, 1), each = ncol(object))
    groupFactor <- factor(groupValues, levels = c(0, 1))

    # Build class object
    exptData <- SummarizedExperiment(
        assays = assaysList,
        rowRanges = rowRanges(object),
        colData = DataFrame(groupFactor = groupFactor)
    )

    # Return class object
    exptData
}


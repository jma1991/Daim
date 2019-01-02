#' Brand-Altman plot
#'
#' Produce a plot of intensity versus difference.
#'
#' @param object A \code{RangedSummarizedExperiment} object.
#' @param contrast A character vector with exactly two elements: the name or index of the numerator level for the fold change, and the name or index of the denominator level for the change.
#' #' @return A Brand-Altman plot on the current graphics device.

plotMA <- function(object, contrast = NULL) {

    # Check argument missing
    if (missing(object)) {
        stop("`object` is missing, with no default.", call. = FALSE)
    }
    if (missing(contrast)) {
        stop("`contrast` is missing, with no default.", call. = FALSE)
    }

    # Check argument class
    if (class(object) != "RangedSummarizedExperiment") {
        stop("`object` must be a RangedSummarizedExperiment class.", call. = FALSE)
    }

    # Define candidate values
    libNames <- c(seq_len(ncol(object)), colnames(object))

    # Check first contrast value
    if (!(contrast[1] %in% libNames)) {
        stop("Can't find column `", contrast[1], "` in `object`.", call. = FALSE)
    }

    # Check second contrast value
    if (!(contrast[2] %in% libNames)) {
        stop("Can't find column `", contrast[2], "` in `object`.", call. = FALSE)
    }

    # Dispatch relevant method
    UseMethod("plotMA", object)

}

plotMA.RangedSummarizedExperiment <- function(object, contrast) {

    # Filter by size and methylation
    sizeFilter <- mcols(object)$filterBySize
    exprFilter <- mcols(object)$filterByExpr

    # Retrieve assay matrix
    if (is.logical(sizeFilter) & is.logical(exprFilter)) {
        assayFilter <- sizeFilter & exprFilter
        assayMatrix <- assay(object[assayFilter])
    } else {
        assayMatrix <- assay(object)
    }

    # Get contrast values
    colIndex1 <- contrast[1]
    colIndex2 <- contrast[2]

    # Calculate mean difference values
    sideValues1 <- rowMeans(assayMatrix[, contrast])
    sideValues2 <- assayMatrix[, colIndex1] - assayMatrix[, colIndex2]

    # Calculate loess model fit
    modelFit <- limma::loessFit(x = sideValues1, y = sideValues2)
    orderFit <- order(sideValues1)

    # Create pretty breaks
    sideBreaks1 <- pretty(sideValues1, n = 10)
    sideBreaks2 <- pretty(sideValues2, n = 10)

    # Create pretty limits
    sideLimits1 <- range(sideBreaks1)
    sideLimits2 <- range(sideBreaks2)

    # Create sample labels
    nullNames <- is.null(colnames(object))
    if (nullNames) {
        sampleNames <- paste("Sample", contrast)
    } else {
        sampleNames <- colnames(object)[contrast]
    }

    # Define title label
    mainText <- paste(sampleNames, collapse = " vs ")

    # Define colour palette
    brewerPalette <- c("#FFFFFF", RColorBrewer::brewer.pal(n = 5, name = "Greys"))
    colourPalette <- colorRampPalette(brewerPalette)

    # Draw plotting area
    smoothScatter(x = sideValues1,
                  y = sideValues2,
                  colramp = colourPalette,
                  main = mainText,
                  nrpoints = 100,
                  xlab = "Mean value (A)",
                  ylab = "Log ratio (M)",
                  xlim = sideLimits1,
                  ylim = sideLimits2,
                  xaxt = "n",
                  yaxt = "n")

    # Draw guide line
    abline(h = 0, lty = 2, col = "#e41a1c")

    # Draw model line
    lines(sideValues1[orderFit], modelFit$fitted[orderFit], col = "#4daf4a")

    # Draw axis ticks
    axis(side = 1, at = sideBreaks1, las = 0)
    axis(side = 2, at = sideBreaks2, las = 2)

}

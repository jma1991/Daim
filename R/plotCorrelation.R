#' Plot correlation between two samples
#'
#' Produce a correlation plot between two samples.
#'
#' @param object A \code{RangedSummarizedExperiment} object.
#' @param contrast The name or index of two samples in the object.

plotCorrelation <- function(object, contrast = NULL) {

    # Check argument missing
    if (missing(object)) {
        stop("`object` is missing, with no default.", call. = FALSE)
    }
    if (missing(contrast)) {
        stop("`contrast` is missing, with no default.", call. = FALSE)
    }

    # Check argument class
    if (!is(object, "RangedSummarizedExperiment")) {
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
    UseMethod("plotCorrelation", object)
}

plotCorrelation.RangedSummarizedExperiment <- function(object, contrast) {

    # Filter by size and methylation
    sizeFilter <- mcols(object)$filterBySize
    exprFilter <- mcols(object)$filterByExpr

    # Retrieve assay matrix
    if (is.logical(sizeFilter) & is.logical(exprFilter)) {
        object <- object[sizeFilter & exprFilter]
    }

    # Get contrast values
    colIndex1 <- contrast[1]
    colIndex2 <- contrast[2]

    # Get assay values
    sideValues1 <- assay(object)[, colIndex1]
    sideValues2 <- assay(object)[, colIndex2]

    # Calculate spearman correlation
    rhoValue <- cor(sideValues1, sideValues2, method = "spearman")
    rhoTitle <- sprintf("Spearman's rho: %.2f", round(rhoValue, digits = 2))

    # Create pretty breaks
    sideBreaks1 <- pretty(sideValues1)
    sideBreaks2 <- pretty(sideValues2)

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

    # Define colour palette
    brewerPalette <- c("#FFFFFF", RColorBrewer::brewer.pal(n = 5, name = "Greys"))
    colourPalette <- grDevices::colorRampPalette(brewerPalette)

    # Draw plotting area
    smoothScatter(
        x = sideValues1,
        y = sideValues2,
        colramp = colourPalette,
        main = rhoTitle,
        nrpoints = 100,
        xlab = sampleNames[1],
        ylab = sampleNames[2],
        xlim = sideLimits1,
        ylim = sideLimits2,
        xaxt = "n",
        yaxt = "n")

    # Draw guide line
    abline(0, 1, lty = 2, col = "#e41a1c")

    # Draw axis ticks
    axis(side = 1, at = sideBreaks1, las = 0)
    axis(side = 2, at = sideBreaks2, las = 2)

}

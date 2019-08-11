#' Plot Abundance Densities
#'
#' Plot the density of abundance values for multiple samples on the same plot.
#'
#' @param object A \code{matrix} or \code{RangedSummarizedExperiment} object.
#' @param group A factor specifying group colours.
#' @return An abundance densitites plot on the current graphics device.

plotDensity <- function(object, group = NULL) {

    # Check argument missing
    if (missing(object)) {
        stop("`object` is missing, with no default.", call. = FALSE)
    }

    # Check argument class
    if (!is(object, "RangedSummarizedExperiment")) {
        stop("`object` must be a RangedSummarizedExperiment class.", call. = FALSE)
    }

    # Dispatch relevant method
    UseMethod("plotDensity", object)
}

plotDensity.RangedSummarizedExperiment <- function(object, group = NULL) {

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

    # Compute kernel density
    densityList <- apply(assayMatrix, 2, density)

    # Define pretty breaks
    sideBreaks1 <- pretty(unlist(vapply(densityList, "[", "x", numeric(1))))
    sideBreaks2 <- pretty(unlist(vapply(densityList, "[", "y", numeric(1))))

    # Define pretty limits
    sideLimits1 <- range(sideBreaks1)
    sideLimits2 <- range(sideBreaks2)

    # Create colour palette
    colorPalette <- RColorBrewer::brewer.pal(n = 9, name = "Set1")

    # Define colours by group
    if (is.null(group)) {
        groupFactor <- as.factor(seq_len(ncol(object)))
        groupColors <- rep(colorPalette, length.out = length(groupFactor))
    } else {
        groupFactor <- as.factor(group)
        groupColors <- colorPalette[as.integer(groupFactor)]
    }

    # Draw plotting area
    plot(
        NA,
        type = "n",
        axes = TRUE,
        main = "Read density",
        xlab = "Abundance",
        ylab = "Density",
        xlim = sideLimits1,
        ylim = sideLimits2,
        xaxt = "n",
        yaxt = "n")

    # Draw axis ticks
    axis(side = 1, at = sideBreaks1, las = 0)
    axis(side = 2, at = sideBreaks2, las = 2)

    # Draw kernel density lines
    mapply(lines, densityList, col = groupColors)

    # Draw plot legend
    legend(
        "topright",
        legend = levels(groupFactor),
        col = unique(groupColors),
        bty = "n",
        lty = 1)

}

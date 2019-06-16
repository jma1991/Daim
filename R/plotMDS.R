#' Multidimensional Scaling
#'
#' Plot samples on principal coordinates.
#'
#' @param object A \code{matrix} or \code{RangedSummarizedExperiment} object.
#' @param group A factor specifying group colours.
#' @return A multidimensional scaling plot on the current graphics device.

plotMDS <- function(object, group) {

    # Check argument missing
    if (missing(object)) {
        stop("`object` is missing, with no default.", call. = FALSE)
    }

    # Check argument class
    if (!(class(object) %in% c("matrix", "RangedSummarizedExperiment"))) {
        stop("`object` must be a matrix or RangedSummarizedExperiment class.", call. = FALSE)
    }

    # Dispatch method
    UseMethod("plotMDS", object)

}

plotMDS.matrix <- function(object, group = NULL) {

    # Retrieve assay matrix
    assayMatrix <- object

    # Remove ranges with zero variance
    keepRowVars <- which(matrixStats::rowVars(assayMatrix) != 0)
    assayMatrix <- assayMatrix[keepRowVars, ]

    # Multidimensional scaling
    assayDists <- dist(t(assayMatrix))
    assayScale <- cmdscale(assayDists, eig = TRUE, k = 2)

    # Define plotting data
    sideValues1 <- assayScale$points[, 1]
    sideValues2 <- assayScale$points[, 2]

    # Define pretty breaks
    sideBreaks1 <- pretty(sideValues1)
    sideBreaks2 <- pretty(sideValues2)

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
    plot(x = sideValues1,
         y = sideValues2,
         pch = 21,
         cex = 2,
         bg = groupColors,
         main = "Classical multidimensional scaling",
         xlab = "Coordinate 1",
         ylab = "Coordinate 2",
         xlim = sideLimits1,
         ylim = sideLimits2,
         xaxt = "n",
         yaxt = "n")

    # Draw axis ticks
    axis(side = 1, at = sideBreaks1, las = 0)
    axis(side = 2, at = sideBreaks2, las = 2)

    # Draw plot legend
    legend("topleft",
           legend = levels(groupFactor),
           pch = 21,
           cex = 1,
           pt.bg = unique(groupColors),
           bty = "n")

}

plotMDS.RangedSummarizedExperiment <- function(object, group = NULL) {

    # Filter by size and methylation
    sizeFilter <- mcols(object)$filterBySize
    exprFilter <- mcols(object)$filterByExpr

    # Retrieve assay matrix
    if (is.null(sizeFilter) | is.null(exprFilter)) {
        assayMatrix <- assay(object, "qsmoothData")
    } else {
        assayFilter <- sizeFilter & exprFilter
        assayMatrix <- assay(object, "qsmoothData")[assayFilter, ]
    }

    # Multidimensional scaling
    assayDists <- dist(t(assayMatrix))
    assayScale <- cmdscale(assayDists, eig = TRUE, k = 2)

    # Define plotting data
    sideValues1 <- assayScale$points[, 1]
    sideValues2 <- assayScale$points[, 2]

    # Define pretty breaks
    sideBreaks1 <- pretty(sideValues1)
    sideBreaks2 <- pretty(sideValues2)

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
    plot(x = sideValues1,
         y = sideValues2,
         pch = 21,
         cex = 2,
         bg = groupColors,
         main = "Multidimensional scaling",
         xlab = "Coordinate 1",
         ylab = "Coordinate 2",
         xlim = sideLimits1,
         ylim = sideLimits2,
         xaxt = "n",
         yaxt = "n")

    # Draw axis ticks
    axis(side = 1, at = sideBreaks1, las = 0)
    axis(side = 2, at = sideBreaks2, las = 2)

    # Draw plot legend
    legend("topleft",
           legend = levels(groupFactor),
           pch = 21,
           cex = 1,
           pt.bg = unique(groupColors),
           bty = "n")

}

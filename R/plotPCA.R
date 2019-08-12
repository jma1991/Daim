#' Principal Components Analysis
#'
#' Performs a principal components analysis on the given data and plots the results.
#'
#' @param object A \code{matrix} or \code{RangedSummarizedExperiment} object.
#' @param group A factor specifying group colours.
#' @return A principal components plot on the current graphics device.

plotPCA <- function(object, group = NULL) {

    # Check argument missing
    if (missing(object)) {
        stop("`object` is missing, with no default.", call. = FALSE)
    }

    # Check argument class
    if (!(class(object) %in% c("matrix", "RangedSummarizedExperiment"))) {
        stop("`object` must be a matrix or RangedSummarizedExperiment class.", call. = FALSE)
    }

    # Dispatch relevant method
    UseMethod("plotPCA", object)

}

plotPCA.matrix <- function(object, group = NULL) {

    # Retrieve assay matrix
    assayMatrix <- object

    # Remove ranges with zero variance
    keepRowVars <- which(matrixStats::rowVars(assayMatrix) != 0)
    assayMatrix <- assayMatrix[keepRowVars, ]

    # Identify principal components
    compData <- prcomp(t(assayMatrix), center = TRUE, scale. = TRUE)
    eigValue <- compData$sdev^2
    compVars <- round((eigValue / sum(eigValue)) * 100, digits = 2)

    # Define plotting data
    sideValues1 <- compData$object[, 1]
    sideValues2 <- compData$object[, 2]

    # Define pretty breaks
    sideBreaks1 <- pretty(sideValues1)
    sideBreaks2 <- pretty(sideValues2)

    # Define pretty limits
    sideLimits1 <- range(sideBreaks1)
    sideLimits2 <- range(sideBreaks2)

    # Create colour palette
    colorPalette <- RColorBrewer::brewer.pal(n = 3, name = "Dark2")

    # Define colours by group
    if (is.null(group)) {
        groupFactor <- as.factor(seq_len(ncol(object)))
        groupColors <- rep(colorPalette, length.out = length(groupFactor))
    } else {
        groupFactor <- as.factor(group)
        groupColors <- colorPalette[as.integer(groupFactor)]
    }

    # Draw plotting area
    plot(object = sideValues1,
         y = sideValues2,
         pch = 21,
         bg = groupColors,
         main = "Principal components",
         xlab = paste0("PC1 (", compVars[1], "%)"),
         ylab = paste0("PC2 (", compVars[2], "%)"),
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
           pt.bg = unique(groupColors),
           bty = "n")
}

plotPCA.RangedSummarizedExperiment <- function(object, group = NULL) {

    # Filter by size and methylation
    sizeFilter <- mcols(object)$filterBySize
    exprFilter <- mcols(object)$filterByExpr

    # Retrieve assay matrix
    if (is.logical(sizeFilter) & is.logical(exprFilter)) {
        assayFilter <- sizeFilter & exprFilter
        assayMatrix <- assay(object, "qsmoothData")[assayFilter, ]
    } else {
        assayMatrix <- assay(object, "qsmoothData")
    }

    # Identify principal components
    compData <- prcomp(t(assayMatrix), center = TRUE, scale. = TRUE)
    eigValue <- compData$sdev^2
    compVars <- round((eigValue / sum(eigValue)) * 100, digits = 2)

    # Define plotting data
    sideValues1 <- compData$x[, 1]
    sideValues2 <- compData$x[, 2]

    # Define pretty breaks
    sideBreaks1 <- pretty(sideValues1)
    sideBreaks2 <- pretty(sideValues2)

    # Define pretty limits
    sideLimits1 <- range(sideBreaks1)
    sideLimits2 <- range(sideBreaks2)

    # Create colour palette
    colorPalette <- RColorBrewer::brewer.pal(n = 3, name = "Dark2")

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
         main = "Principal components",
         xlab = paste0("PC1 (", compVars[1], "%)"),
         ylab = paste0("PC2 (", compVars[2], "%)"),
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

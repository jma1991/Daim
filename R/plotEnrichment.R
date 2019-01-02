#' Distribution of reads along the genome
#'
#' Produce a cumulative sum curve showing the fraction of reads aligned to all restriction fragments in the genome.
#'
#' @param object A \code{RangedSummarizedExperiment} object.
#' @param group A factor specyfing columns of object to groups.
#' @return An enrichment plot on the current graphics device.

plotEnrichment <- function(object, group = NULL) {

    # Check argument missing
    if (missing(object)) {
        stop("`object` is missing, with no default.", call. = FALSE)
    }

    # Check argument class
    if (!is(object, "RangedSummarizedExperiment")) {
        stop("`object` must be a RangedSummarizedExperiment class.", call. = FALSE)
    }

    # Check argument value
    if (!("countsData" %in% assayNames(object))) {
        stop("Can't find matrix 'countsData' in assays slot of 'object'.", call. = FALSE)
    }

    # Dispatch relevant method
    UseMethod("plotEnrichment", object)
}

plotEnrichment.RangedSummarizedExperiment <- function(object, group = NULL) {

    # Downsample counts matrix
    countsData <- assay(object, "countsData")
    minLibSize <- min(colSums(countsData))
    propValues <- minLibSize / colSums(countsData)
    sampleData <- DropletUtils::downsampleMatrix(countsData, prop = propValues)

    # Compute cumulative ranks
    rankValues <- seq_len(nrow(sampleData)) / nrow(sampleData)
    fracValues <- apply(sampleData, 2, function(x) cumsum(sort(x)) / max(cumsum(sort(x))))

    # Downsample for efficiency
    keepValues <- seq(1, nrow(sampleData), length.out = 100)
    rankValues <- rankValues[keepValues]
    fracValues <- fracValues[keepValues, ]

    # Define pretty breaks
    sideBreaks <- seq(0, 1, by = 0.2)

    # Define pretty limits
    sideLimits <- range(sideBreaks)

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
    matplot(x = rankValues,
            y = fracValues,
            type = "l",
            lty = 1,
            col = groupColors,
            main = "Sample enrichment",
            xlab = "Fraction of reads",
            ylab = "Fraction of fragments",
            xlim = sideLimits,
            ylim = sideLimits,
            xaxt = "n",
            yaxt = "n")

    # Draw axis ticks
    axis(side = 1, at = sideBreaks, las = 0)
    axis(side = 2, at = sideBreaks, las = 2)

    # Draw uniform line
    abline(0, 1, lty = 2)

    # Draw plot legend
    legend("topleft",
           legend = levels(groupFactor),
           col = unique(groupColors),
           bty = "n",
           lty = 1)

}

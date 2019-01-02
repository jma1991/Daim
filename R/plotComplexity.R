#' Plot restriction fragment complexity
#'
#' Produce a restriction fragment complexity plot.
#'
#' @param object A RangedSummarizedExperiment object.
#' @param group A grouping factor.
#' @return A restriction fragment complexity plot on the current graphics device.

plotComplexity <- function(object, group = NULL) {

    # Check argument missing
    if (missing(object)) {
        stop("`object` is missing, with no default.", call. = FALSE)
    }

    # Check argument class
    if (class(object) != "RangedSummarizedExperiment") {
        stop("`object` must be a RangedSummarizedExperiment class.", call. = FALSE)
    }

    # Check argument value
    if (!("countsData" %in% assayNames(object))) {
        stop("Can't find matrix 'countsData' in assays slot of 'object'.", call. = FALSE)
    }

    # Dispatch relevant method
    UseMethod("plotComplexity", object)
}

plotComplexity.RangedSummarizedExperiment <- function(object, group = NULL) {

    # Calculate complexity curve for each library
    curveList <- apply(assay(object), 2, function(libCounts) {

        # Compute library size
        libSize <- sum(libCounts)

        # Define step size in proportions
        stepSize <- seq(0, libSize, by = 1e6)

        # Compute proportions vector
        propSeq <- stepSize / libSize

        # Correct proportions vector
        propSeq <- c(propSeq[propSeq < 1], 1)

        # Define proportions length
        propLen <- length(propSeq)

        # Downsample counts matrix
        repMatrix <- matrix(rep(libCounts, each = propLen), ncol = propLen, byrow = TRUE)
        repMatrix <- DropletUtils::downsampleMatrix(repMatrix, prop = propSeq, bycol = TRUE)

        # Compute total reads and distinct fragments
        totalReads <- colSums(repMatrix) / 1e6
        totalFrags <- colSums(repMatrix > 0) / 1e6

        # Return complexity curve data
        data.frame(totalReads, totalFrags)

    })

    # Define max limits
    maxReads <- max(sapply(curveList, function(x) max(x$totalReads, na.rm = TRUE)))
    maxFrags <- max(sapply(curveList, function(x) max(x$totalFrags, na.rm = TRUE)))

    # Create pretty breaks
    sideBreaks1 <- pretty(seq(0, maxReads, by = 1))
    sideBreaks2 <- pretty(seq(0, maxFrags, length.out = 5))

    # Create pretty limits
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
        x = NA,
        y = NA,
        type = "l",
        lty = 1,
        axes = TRUE,
        main = "Fragment complexity",
        xlab = "Total reads (M)",
        ylab = "Distinct fragments (M)",
        xlim = sideLimits1,
        ylim = sideLimits2,
        xaxt = "n",
        yaxt = "n"
    )

    # Draw complexity curves
    mapply(lines, curveList, col = groupColors)

    # Draw axis ticks
    axis(side = 1, at = sideBreaks1, las = 0)
    axis(side = 2, at = sideBreaks2, las = 2)

    # Draw plot legend
    legend(
        x = "bottomright",
        legend = levels(groupFactor),
        col = unique(groupColors),
        bty = "n",
        lty = 1
    )

}

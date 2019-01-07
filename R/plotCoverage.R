#' Cumulative distribution of sequencing coverage
#'
#' Produce a cumulative distribution of sequence coverge plot.
#'
#' @param object A \code{RangedSummarizedExperiment} object.
#' @param group A factor assigning samples to groups.
#' @return A cumulative distribution of sequence coverge plot on the current graphics device.

plotCoverage <- function(object, group = NULL) {

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
    UseMethod("plotCoverage", object)
}

plotCoverage.RangedSummarizedExperiment <- function(object, group = NULL) {

    # Compute cumulative distribution function
    distFun <- apply(assay(object, "countsData"), 2, ecdf)
    distSeq <- seq(0, 100, length.out = 100)
    distCov <- lapply(distFun, function(x) x(distSeq))
    distCov <- 1 - t(do.call(rbind, distCov))

    # Define pretty breaks
    sideBreaks1 <- pretty(distSeq)
    sideBreaks2 <- pretty(seq(0, max(distCov) + 0.1, by = 0.1))

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
    matplot(x = distSeq,
            y = distCov,
            type = "l",
            lty = 1,
            col = groupColors,
            axes = TRUE,
            main = "Read coverage",
            xlab = "Read count",
            ylab = "Fraction of fragments >= Read count",
            xlim = sideLimits1,
            ylim = sideLimits2,
            xaxt = "n",
            yaxt = "n")

    # Draw axis ticks
    axis(side = 1, at = sideBreaks1, las = 0)
    axis(side = 2, at = sideBreaks2, las = 2)

    # Draw plot legend
    legend("topright",
           legend = levels(groupFactor),
           col = unique(groupColors),
           bty = "n",
           lty = 1)
}

#' plotTSS
#'
#' Plot the distance between peak regions and their putatively regulated genes.
#'
#' @param object A \code{GRanges} object of annotated peak regions.

plotTSS <- function(object, mode = c("absolute", "relative")) {

    # Check argument missing
    if (missing(object)) {
        stop("`object` is missing, with no default.", call. = FALSE)
    }

    # Check argument class
    if (!is(object, "GRanges")) {
        stop("`object` must be of class GRanges.", call. = FALSE)
    }

    # Check argument value
    if (is.null(object$distanceTSS)) {
        stop("Column `distanceTSS` must exist in metadata columns of `object`.", call. = FALSE)
    }

    # Check chosen mode is valid
    chosenMode <- match.arg(mode)

    # Dispatch relevant method
    if (chosenMode == "absolute") {
        plotTSS.absolute(object)
    } else if (chosenMode == "relative") {
        plotTSS.relative(object)
    } else {
        stop('`mode` should be one of "absolute" or "relative".', call. = FALSE)
    }

}

plotTSS.absolute <- function(object) {

    # Compute absolute distance
    absDist <- abs(mcols(object)$distanceTSS)

    # Divide distance into intervals
    cutBreaks <- c(0, 5, 50, 500, Inf) * 1e3
    cutRanges <- cut(absDist, cutBreaks, right = FALSE)
    cutCounts <- table(cutRanges)
    cutValues <- (cutCounts / sum(cutCounts)) * 100

    # Create pretty labels
    cutLabels <- c("0 to 5", "5 to 50", "50 to 500", "> 500")
    names(cutValues) <- cutLabels

    # Create colour palette
    colorPalette <- RColorBrewer::brewer.pal(n = 5, name = "Accent")[5]

    # Create accurate Y axis
    axisSpace <- 1.1
    axisLimit <- range(pretty(c(0, cutValues))) * axisSpace

    # Generate plot of binned distances
    myPlot <- barplot(
        height = cutValues,
        col = colorPalette,
        las = 1,
        main = "Binned by absolute distance to TSS",
        xlab = "Absolute distance to TSS (kb)",
        ylab = "Peak-TSS annotations (%)",
        ylim = axisLimit
    )

    # Add frequency above each bar
    text(x = myPlot, y = cutValues, label = cutCounts, pos = 3, cex = 0.75)

}

plotTSS.relative <- function(object) {

    # Compute relative distance
    relDist <- mcols(object)$distanceTSS

    # Divide distance into intervals
    cutBreaks <- c(-Inf, -500, -50, -5, 0, 5, 50, 500, Inf) * 1e3
    cutRanges <- cut(relDist, cutBreaks, right = FALSE)
    cutCounts <- table(cutRanges)
    cutValues <- (cutCounts / sum(cutCounts)) * 100

    # Create pretty labels
    cutLabels <- c("< -500", "-500 to -50", "-50 to -5", "-5 to 0", "0 to 5", "5 to 50", "50 to 500", "> 500")
    names(cutValues) <- cutLabels

    # Create colour palette
    colorPalette <- RColorBrewer::brewer.pal(n = 5, name = "Accent")[5]

    # Create accurate Y axis
    axisSpace <- 1.1
    axisLimit <- range(pretty(c(0, cutValues))) * axisSpace

    # Generate plot of binned distances
    myPlot <- barplot(
        height = cutValues,
        col = colorPalette,
        las = 1,
        main = "Binned by relative distance to TSS",
        xlab = "Relative distance to TSS (kb)",
        ylab = "Peak-TSS annotations (%)",
        ylim = axisLimit
    )

    # Add frequency above each bar
    text(x = myPlot, y = cutValues, label = cutCounts, pos = 3, cex = 0.75)

}

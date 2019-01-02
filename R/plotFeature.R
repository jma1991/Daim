#' plotFeature
#'
#' @param object A \code{GRanges} object of annotated peak regions.
#'

plotFeature <- function(object) {

    # Tabulate annotations
    cutLabels <- mcols(object)$genomicFeature
    cutLevels <- c("TSS", "TES", "Exon", "UTR5", "UTR3", "Intron", "Intergenic")
    cutCounts <- table(factor(cutLabels, levels = cutLevels))
    cutValues <- (cutCounts / sum(cutCounts)) * 100

    # Sort by desired order
    cutCounts <- cutCounts[cutLevels]
    cutValues <- cutValues[cutLevels]

    # Create colour palette
    colorPalette <- RColorBrewer::brewer.pal(n = 7, name = "Accent")

    # Create accurate Y axis
    axisSpace <- 1.1
    axisLimit <- range(pretty(c(0, cutValues))) * axisSpace

    # Generate plot of binned distances
    myPlot <- barplot(
        height = cutValues,
        col = colorPalette,
        las = 1,
        main = "Binned by genomic feature",
        xlab = "Genomic feature",
        ylab = "Peak annotations (%)",
        ylim = axisLimit
    )

    # Add frequency above each bar
    text(x = myPlot, y = cutValues, label = cutCounts, pos = 3, cex = 0.75)

}

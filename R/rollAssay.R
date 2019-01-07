#' Summarize counts from neighbouring ranges
#'
#' Compute the average count from neighbouring ranges.
#'
#' @param object A RangedSummarizedExperiment object.
#' @param size Size allowed between ranges.
#' @return A matrix of summarized counts.

rollAssay <- function(object, size = NULL) {

    # Build container for matrices
    chromLevels <- seqlevels(object)
    chromCount <- length(chromLevels)
    matrixList <- vector("list", chromCount)
    matrixList <- setNames(matrixList, chromLevels)

    # Split assay by chromosome
    chromNames <- seqnames(object)
    assaysList <- split(object, chromNames)

    # Append values to matrices
    for (chromLevel in chromLevels) {

        # Find neighbouring fragments
        queryRanges <- rowRanges(assaysList[[chromLevel]])
        hitsObject <- findOverlaps(queryRanges, queryRanges, maxgap = size)

        # Extract indices
        queryIndex <- queryHits(hitsObject)
        queryCount <- table(queryIndex)
        subjectIndex <- subjectHits(hitsObject)

        # Calculate average count
        assayData <- assay(assaysList[[chromLevel]])
        meanCount <- sweep(rowsum(assayData[subjectIndex, ], queryIndex), 1, queryCount, FUN = "/")

        # Append count to container
        matrixList[[chromLevel]] <- meanCount

    }

    # Return binned average
    baseMatrix <- do.call(rbind, matrixList)

    # Return binned average
    baseMatrix
}

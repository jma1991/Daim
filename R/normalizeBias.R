#' Normalize restriction fragment biases
#'
#' Perform normalization to remove restriction fragment biases.
#'
#' @param object A RangedSummarizedExperiment object.
#' @param group A factor specifying Dam and Dam-fusion groups.
#' @return A \code{\link[SummarizedExperiment]{RangedSummarizedExperiment}} object.

normalizeBias <- function(object, group) {

    # Check argument missing
    if (missing(object)) {
        stop("`object` is missing, with no default.", call. = FALSE)
    }
    if (missing(group)) {
        stop("`group` is missing, with no default.", call. = FALSE)
    }

    # Check argument class
    if (!is(object, "RangedSummarizedExperiment")) {
        stop("`object` must be a RangedSummarizedExperiment object.", call. = FALSE)
    }
    if (!is(group, "factor")) {
        stop("`group` must be a factor.", call. = FALSE)
    }

    # Check argument value
    if (ncol(object) != length(group)) {
        stop("Incompatible lengths between `object` and `group`", call. = FALSE)
    }
    if (!setequal(levels(group), c("0", "1"))) {
        stop("`group` factor must have levels '0' and '1' only.", call. = FALSE)
    }

    # Check accessor missing
    if (!"countsData" %in% assayNames(object)) {
        stop("Matrix `countsData` must exist in assays slot of `object`.", call. = FALSE)
    }
    if (!("digestBySize" %in% names(rowData(object)))) {
        stop("Column `digestBySize` must exist in rowData of `object`.", call. = FALSE)
    }
    if (!("digestByProb" %in% names(rowData(object)))) {
        stop("Column `digestByProb` must exist in rowData of `object`.", call. = FALSE)
    }

    # Check accessor class
    if (!(is.numeric(rowData(object)$digestBySize))) {
        stop("Column `digestBySize` must be a numeric vector.", call. = FALSE)
    }
    if (!(is.numeric(rowData(object)$digestByProb))) {
        stop("Column `digestByProb` must be a numeric vector.", call. = FALSE)
    }

    # Dispatch relevant method
    UseMethod("normalizeBias", object)

}

normalizeBias.RangedSummarizedExperiment <- function(object, group) {

    # Conditional quantile normalisation
    cqNorm <- cqn::cqn(assay(object, "countsData"), mcols(object)$digestByProb, mcols(object)$digestBySize, sqn = FALSE)
    cqData <- 2^(cqNorm$y + cqNorm$offset)

    # Counts from abundance transformation
    ctSums <- colSums(assay(object, "countsData")) / colSums(cqData)
    ctNorm <- round(sweep(cqData, 2, ctSums, FUN = "*"))

    # Add sample groups factor to object
    groupFactor <- as.factor(group)
    colData(object)$groupFactor <- group

    # Normalise to log2 counts per million
    ctSize <- colSums(ctNorm) * edgeR::calcNormFactors(ctNorm, method = "TMM", doWeighting = FALSE)
    ctData <- t(log2(t(ctNorm + 0.5)/(ctSize + 1) * 1e6))

    # Smooth quantile normalisation
    qsNorm <- qsmooth::qsmooth(ctData, groupFactor)
    qsData <- qsmooth::qsmoothData(qsNorm)

    # Store filter genes by fragment length
    keepSize <- filterBySize(rowRanges(object), min = 100, max = 100000)
    rowData(object)$filterBySize <- keepSize

    # Store filter genes by expression level
    keepExpr <- filterByExpr(ctNorm, group = groupFactor, lib.size = ctSize, min.count = 10, min.total.count = 50)
    rowData(object)$filterByExpr <- keepExpr

    # Add qsmooth assay data
    colnames(qsData) <- colnames(object)
    rownames(qsData) <- rownames(object)
    assay(object, "qsmoothData") <- qsData

    # Compute base enrichment
    rowData <- rowData(object)
    rowData$baseMean <- rowMeans(qsData)
    rowData$baseMean0 <- rowMeans(qsData[, groupFactor == 0, drop = FALSE])
    rowData$baseMean1 <- rowMeans(qsData[, groupFactor == 1, drop = FALSE])

    # Compute fold enrichment
    rowData$log2FoldChange <- rowData$baseMean1 - rowData$baseMean0

    # Replace range data
    rowData(object) <- as(rowData, "DataFrame")

    # Return class object
    object

}

#' Call peaks from DamID-seq data
#'
#' Test each restriction fragment for differential methylation and combine into peak regions.
#'
#' @param object A RangedSummarizedExperiment object containing normalised counts.
#' @param alpha Threshold adjusted P value for declaring whether a peak is statistically significant.
#' @param lfc A non-negative value which specifies a log2 fold change threshold.
#' @return A GRanges object of peak regions.

callPeaks <- function(object, alpha = 0.05, lfc = 0) {

    # Check argument missing
    if (missing(object)) {
        stop("`object` is missing, with no default.", call. = FALSE)
    }

    # Check argument class
    if (class(object) != "RangedSummarizedExperiment") {
        stop("`object` must be a RangedSummarizedExperiment object.", call. = FALSE)
    }
    if (!is.numeric(alpha)) {
        stop("`alpha` must be of class numeric.", call. = FALSE)
    }
    if (!is.numeric(lfc)) {
        stop("`lfc` must be of class numeric.", call. = FALSE)
    }

    # Check argument value
    if (!(alpha > 0 & alpha < 1)) {
        stop("`alpha` must be between 0 and 1.", call. = FALSE)
    }
    if (lfc < 0) {
        stop("`lfc` must be a positive value.", call. = FALSE)
    }
    if (!"qsmoothData" %in% assayNames(object)) {
        stop("Matrix `qsmoothData` must exist in assays slot of `object`.", call. = FALSE)
    }

    # Check accessor missing
    if (is.null(object$groupFactor)) {
        stop("Column `groupFactor` must exist in colData of `object`.", call. = FALSE)
    }
    if (!("filterByExpr" %in% names(rowData(object)))) {
        stop("Column `filterByExpr` must exist in rowData of `object`.", call. = FALSE)
    }
    if (!("filterBySize" %in% names(rowData(object)))) {
        stop("Column `filterBySize` must exist in rowData of `object`.", call. = FALSE)
    }
    if (!("baseMean" %in% names(rowData(object)))) {
        stop("Column `baseMean` must exist in rowData of `object`.", call. = FALSE)
    }
    if (!("baseMean0" %in% names(rowData(object)))) {
        stop("Column `baseMean0` must exist in rowData of `object`.", call. = FALSE)
    }
    if (!("baseMean1" %in% names(rowData(object)))) {
        stop("Column `baseMean1` must exist in rowData of `object`.", call. = FALSE)
    }

    # Check accessor class
    if (!(is.logical(rowData(object)$filterByExpr))) {
        stop("Column `filterByExpr` must be a logical vector.", call. = FALSE)
    }
    if (!(is.logical(rowData(object)$filterBySize))) {
        stop("Column `filterBySize` must be a logical vector.", call. = FALSE)
    }
    if (!(is.numeric(rowData(object)$baseMean))) {
        stop("Column `baseMean` must be a numeric vector.", call. = FALSE)
    }
    if (!(is.numeric(rowData(object)$baseMean0))) {
        stop("Column `baseMean0` must be a numeric vector.", call. = FALSE)
    }
    if (!(is.numeric(rowData(object)$baseMean1))) {
        stop("Column `baseMean1` must be a numeric vector.", call. = FALSE)
    }

    # Check accessor value
    if (!setequal(levels(object$groupFactor), c("0", "1"))) {
        stop("Column `groupFactor` must have levels '0' and '1' only.", call. = FALSE)
    }

    # Dispatch relevant method
    UseMethod("callPeaks", object)

}

callPeaks.RangedSummarizedExperiment <- function(object, alpha = 0.05, lfc = 0) {

    # Filter by methylation and size
    exprFilter <- rowData(object)$filterByExpr
    sizeFilter <- rowData(object)$filterBySize
    object <- object[exprFilter & sizeFilter]

    # Create design matrix
    designMatrix <- model.matrix(~ 0 + object$groupFactor)
    colnames(designMatrix) <- c("groupFactor0", "groupFactor1")

    # Estimate mean-variance relationship
    qsData <- assay(object, "qsmoothData")
    modelVar <- limma::voomaByGroup(qsData, group = object$groupFactor)

    # Fit linear model for each fragment
    modelFit <- limma::lmFit(modelVar, design = designMatrix)

    # Extract required fragment statistics
    keepData <- c("baseMean", "baseMean0", "baseMean1", "log2FoldChange")
    fragmentTable <- as.data.frame(rowData(object)[, keepData])

    # Compute contrasts from linear model fit
    contrastMatrix <- limma::makeContrasts(groupFactor = "groupFactor1-groupFactor0", levels = designMatrix)
    modelFit <- limma::contrasts.fit(modelFit, contrasts = contrastMatrix)

    # Compute statisitics for differential methylation
    modelFit <- daimTreat(modelFit, lfc = lfc, trend = TRUE, robust = TRUE)

    # Extract tests from linear model fit
    fragmentTests <- limma::topTable(modelFit, coef = "groupFactor", number = Inf, sort.by = "none")

    # Filter fragments by fold change
    filterByTests <- fragmentTests$logFC >= 0
    fragmentTests <- fragmentTests[filterByTests, ]
    fragmentTable <- fragmentTable[filterByTests, ]

    # Combine fragments into windows
    fragmentRanges <- rowRanges(object)[filterByTests]
    windowRanges <- csaw::mergeWindows(fragmentRanges, tol = 100)

    # Combine statistics across multiple tests
    windowTests <- csaw::combineTests(windowRanges$id, fragmentTests, pval.col = "P.Value", fc.col = "logFC")

    # Calculate statistics for each window
    windowTable <- split(fragmentTable, windowRanges$id)
    windowTable <- lapply(windowTable, colMeans)
    windowTable <- as.data.frame(do.call(rbind, windowTable))

    # Adjust significance by IHW
    windowTable$pval <- windowTests$P.Value
    ihwResult <- IHW::ihw(windowTable$pval, windowTable$baseMean, alpha = alpha)
    windowTable$padj <- IHW::adj_pvalues(ihwResult)

    # Create window ranges
    windowRanges <- windowRanges$region
    mcols(windowRanges) <- windowTable

    # Filter window ranges by alpha
    filterByAlpha <- windowRanges$padj < alpha
    windowRanges <- windowRanges[filterByAlpha, ]
    seqinfo(windowRanges) <- seqinfo(object)

    # Return ranges object
    windowRanges

}

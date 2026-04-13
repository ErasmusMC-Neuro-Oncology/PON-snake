# edited poolRuns for compatibility
poolRuns <- function(object, samples, force=FALSE) {
    library(matrixStats)
    if (!force && "segmented" %in% assayDataElementNames(object))
        stop("Data has already been segmented. Pooling runs will ",
            "remove segmentation (and possible calling) ",
            "results. Please specify force=TRUE, if you want this.")

    if (length(samples) != ncol(object))
        stop("Parameter samples must be a vector of equal length as there are ",
            " samples in object.")

    newsamples <- sort(unique(samples))
    if (length(newsamples) == ncol(object)) {
        sampleNames(object) <- samples
        object <- object[,order(samples)]
        return(object)
    }
    oldsamples <- sampleNames(object)

    bins <- featureData(object)
    phenodata <- pData(object)
    phenodata[, 1] <- samples
    newphenodata <- phenodata[0, ]
    if (inherits(object, "QDNAseqReadCounts")) {
        counts <- assayDataElement(object, "counts")
        newcounts <- matrix(NA_integer_,
            nrow=nrow(object), ncol=length(newsamples),
            dimnames=list(featureNames(object), newsamples))
    } else if (inherits(object, "QDNAseqCopyNumbers")) {
        copynumber <- assayDataElement(object, "copynumber")
        newcopynumber <- matrix(NA_real_,
            nrow=nrow(object), ncol=length(newsamples),
            dimnames=list(featureNames(object), newsamples))
    }

    concatenateIfNotEqual <- function(x) {
        x <- sort(unique(x))
        paste(x, collapse=";")
    }

    for (newsample in newsamples) {
        replicates <- samples == newsample
        if (inherits(object, "QDNAseqReadCounts")) {
            newcounts[, newsample] <- rowSums2(counts, cols=replicates)
        } else if (inherits(object, "QDNAseqCopyNumbers")) {
            newcopynumber[, newsample] <- rowMeans2(copynumber, cols=replicates)
        }
        oldphenodata <- phenodata[replicates, ]
        if ("paired.ends" %in% colnames(oldphenodata)) {
            pairedEnds <- unique(oldphenodata$paired.ends)
            if (length(pairedEnds) != 1)
                pairedEnds <- NA_real_
        }
        totalReads <- sum(oldphenodata$total.reads)
        usedReads <- sum(oldphenodata$used.reads)
        numericCols <- sapply(oldphenodata, FUN=is.numeric)
        oldphenodata[1, numericCols] <- colMeans2(as.matrix(oldphenodata[,numericCols]))
        oldphenodata[1, !numericCols] <- apply(oldphenodata[, !numericCols,
            drop=FALSE], MARGIN=2L, concatenateIfNotEqual)
        if ("paired.ends" %in% colnames(oldphenodata))
            oldphenodata[1, "paired.ends"] <- pairedEnds
        oldphenodata[1, "total.reads"] <- totalReads
        oldphenodata[1, "used.reads"] <- usedReads
        newphenodata <- rbind(newphenodata, oldphenodata[1, ])
    }
    rownames(newphenodata) <- "ReferencePool"
    #newphenodata <- newphenodata[,c(1,2,3,4)]
    newphenodata <- AnnotatedDataFrame(newphenodata,
        varMetadata=varMetadata(object))
    if (inherits(object, "QDNAseqReadCounts")) {
        storage.mode(newcounts) <- "integer"
        object2 <- new("QDNAseqReadCounts",
            bins=bins,
            counts=newcounts,
            phenodata=newphenodata)
    } else if (inherits(object, "QDNAseqCopyNumbers")) {
        object2 <- new("QDNAseqCopyNumbers",
            bins=bins,
            copynumber=newcopynumber,
            phenodata=newphenodata)
    }
    expectedVariance <- function(object) {
        binsToUse <- function(object) {
    
    if ("use" %in% colnames(object)) {
        return(object$use)
    } else if ("filter" %in% colnames(object)) {
        object$use <- object$filter
        object$filter <- NULL
        return(object$use)
    }
    rep(TRUE, times=nrow(object))
}
    expected.variance <- sum(binsToUse(object)) / object$used.reads
    if ("paired.ends" %in% colnames(pData(object))) {
        multiplier <- ifelse(object$paired.ends, 2, 1)
        expected.variance <- expected.variance * multiplier
    }
    expected.variance
}
    object2$expected.variance <- expectedVariance(object2)
    object2
}

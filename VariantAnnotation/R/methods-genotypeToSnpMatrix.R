### =========================================================================
### genotypeToSnpMatrix methods 
### =========================================================================

## Coding for snpMatrix :
## 0 = missing OR multiallelic OR multi-ALT values
## 1 = homozygous reference (0|0) 
## 2 = heterozygous (0|1 or 1|0) 
## 3 = homozygous alternate (risk) allele (1|1)

setMethod("genotypeToSnpMatrix", "CollapsedVCF",
          function(x, uncertain=FALSE, ...)
{
    ok <- suppressWarnings(require("snpStats", quietly=TRUE, 
                                   character.only=TRUE))
    ok || stop("'snpStats' required; try biocLite('snpStats')", call.=FALSE) 

    alt <- alt(x)
    if (is(alt, "CompressedCharacterList")) {
        warning("structural variants detected and not supported by SnpMatrix; returning NULL")
        return(NULL)
    }
    ref <- ref(x)

    if (!uncertain) {
        gt <- geno(x)$GT
    } else {
        geno.cols <- row.names(geno(exptData(x)[["header"]]))
        if ("GP" %in% geno.cols) {
            gt <- geno(x)$GP
        } else if ("GL" %in% geno.cols) {
            gt <- GLtoGP(geno(x)$GL)
        } else {
            warning("uncertain=TRUE requires GP or GL; returning NULL")
            return(NULL)
        }
        # it's possible to have a single GL field for ref hom calls
        np <- elementLengths(gp)
        gt[np == 1] <- list(1,0,0)
    }

    callGeneric(gt, ref, alt)
})

setMethod("genotypeToSnpMatrix", "matrix",
          function(x, ref, alt, ...)
{
    # if mode of matrix is character, we have GT with a single value for each snp

    # if mode of matrix is list, we have GP with multiple values for each snp
    # for each sample, call probabilityToSnpMatrix
    # use cbind2 to combine
  
    map <- setNames(sapply(rep(c(0, 1, 2, 2, 3), 2), as.raw),
                    c(".|.", "0|0", "0|1", "1|0", "1|1",
                      "./.", "0/0", "0/1", "1/0", "1/1"))
    diploid <- x %in% names(map)
    if (!all(diploid)) {
        warning("non-diploid variants are set to NA")
        x[!diploid] <- ".|."
    }

    altelt <- elementLengths(alt) == 1L
    if (!all(altelt)) {
        warning("variants with >1 ALT allele are set to NA")
        x[!altelt] <- ".|."
    }

    altseq <- logical(length(alt))
    idx <- rep(altelt, elementLengths(alt))
    altseq[altelt] = width(unlist(alt))[idx] == 1L
    snv <- altseq & (width(ref) == 1L)
    if (!all(snv)) {
        warning("non-single nucleotide variations are set to NA")
        x[!snv] <- ".|."
    }

    mat <- matrix(map[x], nrow=ncol(x), ncol=nrow(x),
                  byrow=TRUE, dimnames=rev(dimnames(x)))
    genotypes <- new("SnpMatrix", mat)

    flt <- !(snv & altelt)
    map <- .createMap(rownames(x), ref, alt, flt)

    list(genotypes = genotypes, map = map)
})

.createMap <- function(nms, ref, alt, flt)
{
    if (is.null(ref))
        DataFrame(snp.names=character(0), 
                  allele.1=DNAStringSet(), 
                  allele.2=DNAStringSetList(),
                  ignore=logical())
    else 
        DataFrame(snp.names=nms, 
                  allele.1=ref, 
                  allele.2=alt,
                  ignore=flt)
}

probabilityToSnpMatrix <- function(probs) {
    if (ncol(probs) != 3)
        stop("input matrix should have 3 columns: P(RR), P(RA), P(AA)")

    # skip missing values when checking for validity of probabilities
    missing <- rowSums(is.na(probs)) > 0
    if (!isTRUE(all.equal(rowSums(probs[!missing,]), rep(1,sum(!missing)),
                          tolerance=0.0001, check.attributes=FALSE,
                          check.names=FALSE)))
        stop("sum of probabilities in each row of input matrix should = 1")

    # post2g can't handle missing data
    if (sum(missing) > 0) {
        probs[missing,] <- 0
        g <- post2g(probs)
        g[missing] <- as.raw(0)
    } else {
        g <- post2g(probs)
    }
    g <-  matrix(g, nrow=1, dimnames=list(NULL,rownames(probs)))
    new("SnpMatrix", g)
}

GLtoGP <- function(gl) {
    gp <- gl
    # there is probably a faster way to do this
    for (i in 1:length(gl)) {
        gp[[i]] <- 10^gl[[i]] / sum(10^gl[[i]])
    }
    gp
}

.matrixOfListsToArray <- function(x) {
    # find number of elements of each cell of x
    n <- elementLengths(x)
    maxn <- max(n)
    
    # for cells with less than the max number of elements, add NAs
    idx <- n < maxn
    x[idx] <- lapply(x[idx], function(a){c(a, rep(NA, maxn-length(a)))})
    
    # unlist and convert to array
    x <- array(unlist(x), dim=c(maxn, nrow(x), ncol(x)),
               dimnames=c(NULL, rownames(x), colnames(x)))
    x <- aperm(x, c(2,3,1))

    x
}

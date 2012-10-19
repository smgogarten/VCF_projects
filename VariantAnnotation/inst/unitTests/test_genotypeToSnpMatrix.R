
test_gSM_matrix_GT <- function() {
    mat <- matrix(c(".|.", "0|0", "0|1", "1|0", "1|1",
                    "./.", "0/0", "0/1", "1/0", "1/1"),
                  ncol=2, dimnames=list(1:5,1:2))
    sm <- new("SnpMatrix",
              matrix(as.raw(c(0, 1, 2, 2, 3,
                              0, 1, 2, 2, 3)),
                     nrow=2, byrow=TRUE, dimnames=list(1:2,1:5)))
    ref <- DNAStringSet(rep("A",5))
    alt <- DNAStringSetList(DNAStringSet("C"),
                            DNAStringSet("G"),
                            DNAStringSet("T"),
                            DNAStringSet("C"),
                            DNAStringSet("G"))
    map <- DataFrame(snp.names=rownames(mat),
                             allele.1=ref, 
                             allele.2=alt,
                             ignore=rep(FALSE,5))
    gtsm <- genotypeToSnpMatrix(mat, ref, alt)
    checkIdentical(sm, gtsm$genotypes)
    checkIdentical(map, gtsm$map)
                             
}

test_gSM_matrix_GT_nondiploid<- function() {
    mat <- matrix(c("0|1", "1|0", "1|1",
                    "1/2", "2/1", "2/2"),
                  ncol=2, dimnames=list(1:3,1:2))
    sm <- new("SnpMatrix",
              matrix(as.raw(rep(0,6)),
                     nrow=2, byrow=TRUE, dimnames=list(1:2,1:3)))
    ref <- DNAStringSet(rep("A",3))
    alt <- DNAStringSetList(DNAStringSet(c("C","G")),
                            DNAStringSet(c("G","T")),
                            DNAStringSet(c("T","C")))
    map <- DataFrame(snp.names=rownames(mat),
                             allele.1=ref, 
                             allele.2=alt,
                             ignore=rep(TRUE,3))
    gtsm <- genotypeToSnpMatrix(mat, ref, alt)
    checkIdentical(sm, gtsm$genotypes)
    checkIdentical(map, gtsm$map)
}

test_gSM_matrix_GT_nonsnv <- function() {
    mat <- matrix(c("0|0", "0|1", "1|0",
                    "0/0", "0/1", "1/0"),
                  ncol=2, dimnames=list(1:3,1:2))
    sm <- new("SnpMatrix",
              matrix(as.raw(rep(0,6)),
                     nrow=2, byrow=TRUE, dimnames=list(1:2,1:3)))
    ref <- DNAStringSet(c("A","ACG","ACG"))
    alt <- DNAStringSetList(DNAStringSet("CGT"),
                            DNAStringSet("G"),
                            DNAStringSet("GAC"))
    map <- DataFrame(snp.names=rownames(mat),
                             allele.1=ref, 
                             allele.2=alt,
                             ignore=rep(TRUE,3))
    gtsm <- genotypeToSnpMatrix(mat, ref, alt)
    checkIdentical(sm, gtsm$genotypes)
    checkIdentical(map, gtsm$map)
}

test_gSM_VCF_structural <- function() {
  fl <- system.file("extdata", "structural.vcf", package="VariantAnnotation")
  vcf <- readVcf(fl, "hg19")
  checkIdentical(NULL, genotypeToSnpMatrix(vcf))
}

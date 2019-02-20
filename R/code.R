#' go2GRangesList
#'
#' Create a GRangesList
#' holding TSS +/- a window of variable length.
#'
#' @export
go2GRangesList <- function(txdb = "TxDb.Hsapiens.UCSC.hg19.knownGene",
                           searchDist = 0,
                           ontologyFilter = c("BP","MF","CC"),
                           evidenceFilter = c("TAS","IC")){
  src <- src_organism(txdb)
  OntFilt <- OntologyFilter(ontologyFilter)
  EvFilt <- EvidenceFilter(evidenceFilter)
  promoters_tbl(src,
                upstream = searchDist, downstream = searchDist,
                columns = "go",
                filter = ~ OntFilt & EvFilt) %>%
    collect %>%
    GRanges %>%
    GenomicRanges::split(~ go) %>%
    return
}


#' matchAnnotations
#'
#' From a RangedSummarizedExperiment, create a new
#' RangedSummarizedExperiment with a sparseMatrix in the assay
#' slot.
#'
#' @export
matchAnnotations <- function(annotations, rse) {
  hits <- findOverlaps(rse,annotations)

  mat <- Matrix::sparseMatrix(i=queryHits(hits),
                             j=subjectHits(hits),
                             x=T,
                             dims=c(queryLength(hits),subjectLength(hits)),
                             dimnames=list(as.character(rowRanges(rse)),names(annotations)))
  stopifnot(rownames(mat)==names(rse))
  rse2 <- SummarizedExperiment(rowRanges(rse),assays = list(annotationMatches=mat %>% set_rownames(NULL)))
  return(rse2)
}

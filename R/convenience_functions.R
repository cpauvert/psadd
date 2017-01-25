#' OTU renaming: from hash to rank abundance prefix
#'
#' OTU names when generated with vsearch for example corresponds
#' to sequence SHA1 hash value. Despite being unique, this notation
#' is not human readable.
#' This function transform these hashes into easier OTU names.
#' Ex: `e443b18` to `OTU_1` for the most abundant OTU
#'
#' This renaming also has the advantage to ease subsetting with OTU names
#' because OTU names are the \code{\link[base]{row.names}} of several
#' \code{\link[phyloseq]{phyloseq}} object.
#'
#' @param physeq \code{\link{phyloseq-class}}
#' @param physeq-renamed \code{\link{phyloseq-class}} with OTU renamed.
#'
#' @importFrom phyloseq taxa_sums taxa_names
#' @export
#'
rename_otu<-function(physeq){
  # OTU Abundance vector
  otu_abd<-taxa_sums(physeq)
  # Get OTU rank after sorting by decreasing abundance
  rank_otus<-names(sort(otu_abd, decreasing = T))
  # Fetch hash otus names
  hash_otus<-taxa_names(physeq)
  # Associate hash to rank abundance prefix
  taxa_names(physeq)<-paste("OTU", match(hash_otus, rank_otus), sep = "_")
  return(physeq)
}

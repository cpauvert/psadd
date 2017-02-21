#' Performs rarefying multiple times and returns averaged OTU table
#'
#' @inheritParams phyloseq::rarefy_even_depth
#' @param n The number of times to perform rarefying of OTU table
#' @param start.seed Random seed to start, it will be incremented up to \code{start.seed+n}
#'
#' @import phyloseq
#'
#' @return An object of class \code{\link[phyloseq]{phyloseq}}.
#'  Only the \code{\link[phyloseq]{phyloseq}} component is modified.
#' @export
#'
#' @seealso \code{\link[phyloseq]{rarefy_even_depth}}
multiple_rrarefy<-function(physeq, sample.size = min(sample_sums(physeq)), n, start.seed){
  newobj<-physeq
  # Initialize OTU table with 0
  otu_table(newobj)<-otu_table(newobj)*0
  for(i in seq(n)){
    otu_table(newobj)<-otu_table(
      otu_table(newobj) + (otu_table(
      rarefy_even_depth(physeq, sample.size,rngseed = start.seed+i,
                        replace = FALSE, trimOTUs = FALSE, verbose = FALSE)) / n),
      taxa_are_rows = taxa_are_rows(physeq))
    message(i," rarefied OTU table out of ",n)
  }
  return(newobj)
}

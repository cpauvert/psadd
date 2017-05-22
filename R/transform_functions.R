#' Performs rarefying multiple times and returns averaged OTU table
#'
#' @inheritParams phyloseq::rarefy_even_depth
#' @param n The number of times to perform rarefying of OTU table
#' @param start.seed Random seed to start, it will be incremented up to \code{start.seed+n}
#'
#' @import phyloseq
#'
#' @return An object of class \code{\link[phyloseq]{phyloseq}}.
#'  Only the \code{\link[phyloseq]{otu_table-class}} component is modified.
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
#' Relevel taxa/OTU to 0 according to a taxa/OTU wise threshold.
#'
#' @param physeq \code{\link[phyloseq]{otu_table-class}}, or \code{\link[phyloseq]{phyloseq-class}}
#' @param otuThreshold A \code{\link{numeric-class} of length equal to the number of taxa/OTU.
#' If the latter is zero, the corresponding taxa/OTU will not be relevel.}
#'
#' @return An object of class \code{\link[phyloseq]{phyloseq}}.
#'  Only the \code{\link[phyloseq]{otu_table-class}} component is modified.
#' @export
#' @references This procedure is described in \href{https://doi.org/10.1128/mSystems.00032-16}{Galan et al. (2016)}.
#'
#'
#' @examples
adaptative_threshold_otu<-function(physeq, otuThreshold){
  # Fetch OTU table
  otu<-otu_table(physeq)
  # set margin depending on OTU orientation: row or column
  taxMargin<-ifelse( taxa_are_rows(otu), 1, 2)
  # remove the number of sequences indicated in otuThreshold
  # by samples for each OTU
  otu.sweeped<-sweep(otu, taxMargin, otuThreshold+1, "-")
  # any OTU abundance negative was hence below the threshold
  # and therefore set to NA
  otu.sweeped[ which(otu.sweeped < 0) ]<-NA
  # Re balance the altered abundance by adding the previously
  # removed threshold. Note that NA are not affected
  otu.sweeped<-sweep(otu.sweeped, taxMargin, otuThreshold+1,"+")
  # OTU that are NA are set to 0 and considered absent.
  otu.sweeped[which(is.na(otu.sweeped))]<-0
  # Return the corrected OTU table
  otu_table(physeq)<-otu_table(otu.sweeped,
                               taxa_are_rows = taxa_are_rows(otu))
  return(physeq)
}

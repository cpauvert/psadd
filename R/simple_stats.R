
#####################################################
#' Apply statistics function to taxa abundance vector
#'
#' It generalize the \code{\link{taxa_sums}} function and enable the use of several summary functions.
#' @param physeq \code{\link{otu_table-class}}, or \code{\link{phyloseq-class}}
#' @param ... The function to be applied to the OTU table
#' @returns A \code{\link{numeric-class}} with length equal to the number of species
#'  in the table, name indicated the taxa ID, and value equal to the result of
#'  the function applied through \code{...}.
#' @seealso \code{\link{taxa_sums}} on which this function is based.
#' @example
#' data(esophagus)
#' taxa_stats(esophagus, mean)[1:6]
#' taxa_stats(esophagus, var)[1:6]
taxa_stats<-function (physeq,...){
  m <- otu_table(physeq)
  marg<-ifelse(taxa_are_rows(m), yes = 1, no = 2)
  return(apply(m, marg, ...))
}

########################################
#' OTU/Species prevalence in all samples
#'
#' Return a named vector of taxa/OTU with the number of samples where
#' their abundance was not null.
#'
#' @param physeq \code{\link{otu_table-class}}, or \code{\link{phyloseq-class}}
#'
#' @return A \code{\link{numeric-class}} with length equal to the number of species
#'  in the table, name indicated the taxa ID, and value equal to the species prevalence
#'  in every samples where it occurs.
#'
#' @seealso \code{\link{taxa_sums}} on which this function is based.
#'  \code{\link{sample_prev}} for a similar output on samples.
#'
#' @examples
#' data(esophagus)
#' taxa_prev(esophagus)
taxa_prev<-function(physeq){
  m <- otu_table(physeq)
  if (taxa_are_rows(m)) {
    rowSums(m > 0)
  }
  else {
    colSums(m > 0)
  }
}


########################################
#' Taxa/OTU number in all samples
#'
#' Return a named vector of samples with the number of taxa/OTU where
#' their abundance was not null. It is similar to the observed richness.
#'
#' @param physeq \code{\link{otu_table-class}}, or \code{\link{phyloseq-class}}
#'
#' @return A \code{\link{numeric-class}} with length equal to the number of samples
#'  in the table, name indicated the sample ID, and value equal to the number of taxa/OTU
#'  each samples contains.
#'
#' @seealso \code{\link{taxa_sums}} on which this function is based.
#'  \code{\link{taxa_prev}} for taxa/OTU prevalence across samples.
#'
#' @examples
#' data(esophagus)
#' sample_prev(esophagus)
sample_prev<-function (physeq){
  m <- otu_table(physeq)
  if (taxa_are_rows(m)) {
    colSums(m > 0)
  }
  else {
    rowSums(m > 0)
  }
}

#########################################
#' Subset samples and prune taxa/OTU with zero abundance
#'
#' @param physeq \code{\link{otu_table-class}}, or \code{\link{phyloseq-class}}
#' @param ... The subsetting expression that should be applied to the
#'  \code{sampleTable}. This is passed on to \code{\link{subset_samples}}.
#'  Note that taxa/OTU with zero abundance are pruned from the \code{\link{phyloseq-class}} object.
#'
#'
#' @return A subsetted object with the same class as \code{physeq}
#' @export
#'
#' @seealso \code{\link{subset_samples}} on which it rely.
#' @examples
#' data(GlobalPatterns)
#' gp<-subset_samples(GlobalPatterns, SampleType=="Ocean")
#' gp.nz<-subset_samples_no_zero(GlobalPatterns, SampleType=="Ocean")
#' ntaxa(gp) # 19216 OTU
#' ntaxa(gp.nz) # 5669 OTU
subset_samples_no_zero<-function(physeq,...){
  # Original subsetting
  x<-subset_samples(physeq,...)
  x<-prune_taxa(
    taxa_sums(x)!=0 , # TRUE/FALSE OTU Vector to prune
    x )
  return(x)
}

#########################################################################
#' Compute taxa/OTU relative abundance
#'
#' @param physeq \code{\link{otu_table-class}}, or \code{\link{phyloseq-class}}
#'
#' @return A \code{\link{numeric-class}} with length equal to the number of species
#'  in the table, name indicated the taxa ID, and value equal to the taxa/OTU
#'  relative abundance
#' @export
#'
#' @examples
relative_abundance<-function(physeq){
  # return relative abundance of OTU in phyloseq object
  abd<-taxa_sums(physeq)
  return(abd/sum(abd))
}

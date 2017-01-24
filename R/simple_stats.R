
#####################################################
#' Apply statistics function to taxa abundance vector
#'
#' It generalize the \code{\link[phyloseq]{taxa_sums}} function and enable the use of several summary functions.
#'
#' @param physeq \code{\link{otu_table-class}}, or \code{\link{phyloseq-class}}
#' @param ... The function to be applied to the OTU table
#' @return A \code{\link{numeric-class}} with length equal to the number of species
#'  in the table, name indicated the taxa ID, and value equal to the result of
#'  the function applied through the passed expression (\code{...}).
#' @seealso \code{\link{taxa_sums}} on which this function is based.
#'
#' @import phyloseq
#'
#' @export
#' @examples
#' require(phyloseq)
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
#' @import phyloseq
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
#####################################################
#' Apply statistics function to samples.
#'
#' It generalize the \code{\link[phyloseq]{sample_sums}} function and enable the use of several summary functions.
#'
#' @param physeq \code{\link{otu_table-class}}, or \code{\link{phyloseq-class}}
#' @param ... The function to be applied to the OTU table
#' @return A \code{\link{numeric-class}} with length equal to the number of samples
#'  in the table. Names indicate sample ID, and values equal to the result of
#'  the function applied through the passed expression (\code{...}).
#' @seealso \code{\link{sample_sums}} on which this function is based.
#'
#' @import phyloseq
#'
#' @export
#' @examples
#' require(phyloseq)
#' data(esophagus)
#' sample_stats(esophagus, mean)[1:6]
#' sample_stats(esophagus, var)[1:6]
sample_stats<-function (physeq,...){
  m <- otu_table(physeq)
  marg<-ifelse(taxa_are_rows(m), yes = 2, no = 1)
  return(apply(m, marg, ...))
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
#' @import phyloseq
#' @examples
#' data(esophagus)
#' sample_prev(esophagus)
sample_prev<-function(physeq){
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
#'
#' @import phyloseq
#'
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
#' Compute taxa/OTU abundance relative to the total dataset
#'
#' Convenience function to divide \code{taxa_sums} output vector by the
#' total number of sequences.
#'
#' @param physeq \code{\link{otu_table-class}}, or \code{\link{phyloseq-class}}
#'
#' @return A \code{\link{numeric-class}} with length equal to the number of species
#'  in the table, name indicated the taxa ID, and value equal to the taxa/OTU
#'  relative abundance
#' @export
#' @import phyloseq
#'
#' @seealso \code{\link{taxa_sums}} on which this function is based.
#' See \code{\link{ctrl_samples_stats}} for a usage example.

relative_abundance<-function(physeq){
  # return relative abundance of OTU in phyloseq object
  abd<-taxa_sums(physeq)
  return(abd/sum(abd))
}

#########################################################################
#' Compute taxa/OTU abundance relative to the total dataset
#'
#' Convenience function to divide \code{taxa_sums} output vector by the
#' total number of sequences.
#' Samples present in the second \code{phyloseq} object --\code{specific}-- are discarded
#'  from the first one --\code{physeq}.
#' Taxa/OTU present in the second \code{phyloseq} object --\code{specific}-- and supposedly
#' of interest are stored.
#' The relative abundances of these taxa/OTU of interest are then computed.
#'
#' @param physeq \code{\link{otu_table-class}}, or \code{\link{phyloseq-class}}
#' @param specific \code{\link{otu_table-class}}, or \code{\link{phyloseq-class}}
#'  containing samples of interest to examine (e.g. Negative or Positive Controls).

#'
#' @return A \code{\link{numeric-class}} with length equal to the number of species
#'  present in the second --\code{specific}-- \code{phyloseq} object.
#'   Names indicate the taxa ID, and value equal to the taxa/OTU
#'  relative abundance in the first --\code{physeq}-- object (pruned of specific samples).
#' @export
#' @import phyloseq
#'
#' @seealso \code{\link{taxa_sums}} on which this function is based.
#' See \code{\link{ctrl_samples_stats}} for a usage example.
relative_abd_without_specific_sample<-function(physeq,specific){
  # Discard sample of interest from the global otu table
  m<-prune_samples(! sample_names(physeq) %in% sample_names(specific),physeq)
  # Keep only OTU in the specific phyloseq object
  otu_to_keep<-taxa_names(specific)
  # Compute Relative Abundance of these OTUs in the first dataset (pruned of specific samples)
  rel_abd<-relative_abundance(m)
  return(rel_abd[ otu_to_keep ])
}

#####################################################################################

#' Generate taxa/OTU statistics to compare control samples with the remaining dataset
#'
#'
#' @inheritParams taxa_stats
#' @param physeq_ctrl \code{\link{otu_table-class}}, or \code{\link{phyloseq-class}}
#'  of control sample. It can be the output of \code{\link{subset_samples_no_zero}}
#'   for example.
#'
#' @details
#' The taxa/OTU prevalence is defined as the number of samples in which it occurs.
#' The prevalence is corrected by substracting the number of samples in the \code{physeq_ctrl}
#' object (usually 1).
#'
#' @return A \code{\link{data.frame}} with the following components
#' \describe{
#' \item{\code{RA_In_Ctrl}}{Taxa/OTU relative abundance in Control sample}
#' \item{\code{RA_In_All}}{Same Taxa/OTU relative abundance in all samples \strong{but} Control samples}
#' \item{\code{Ratio}}{Ratio of \code{RA_In_Ctrl}/\code{RA_In_All}}
#' \item{\code{Prevalence}}{Taxa/OTU Prevalence in all samples \strong{but} Control samples.
#'  Range from 0 to the number of samples in \code{physeq} object
#'   (obtained by \code{\link[phyloseq]{ntaxa}})}
#' }
#'
#' @export
#' @import phyloseq
#' @examples
#' data(soilrep)
#' # Subset a sample
#' soilrep.6CC<-subset_samples_no_zero(soilrep, Sample=="6CC")
#' # Check
#' ctrl_df<-ctrl_samples_stats(soilrep,soilrep.6CC)
#' head(ctrl_df)
ctrl_samples_stats<-function(physeq,physeq_ctrl){
  # RA stands for relative abundance
  RA_ctrl_otu<-relative_abundance(physeq_ctrl)
  RA_without_ctrl_sample<-relative_abd_without_specific_sample(physeq, physeq_ctrl)
  # Prevalence of control OTU in the global dataset, corrected by the number of control sample (1 usually)
  otu_prevalence<-taxa_prev(physeq)[taxa_names(physeq_ctrl)] - nsamples(physeq_ctrl)
  df<-data.frame(
    RA_In_Ctrl=RA_ctrl_otu,
    RA_In_All=RA_without_ctrl_sample,
    Ratio=RA_ctrl_otu/RA_without_ctrl_sample,
    Prevalence=otu_prevalence)
  return(df)
}

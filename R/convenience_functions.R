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
#' @param hash2id_file A character giving the filename to write the
#' correspondence table between the hash and the new ID or "" for output to
#' the console.
#' @return physeq \code{\link{phyloseq-class}} with OTU renamed.
#'
#' @importFrom phyloseq taxa_sums taxa_names
#' @export
#' @examples
#' require(phyloseq)
#' data(esophagus)
#' taxa_names(esophagus)[1:6]
#' taxa_names( rename_otu( esophagus, hash2id_file = "" ) )[1:6]
rename_otu<-function(physeq, hash2id_file){
  # OTU Abundance vector
  otu_abd<-taxa_sums(physeq)
  # Get OTU rank after sorting by decreasing abundance
  rank_otus<-names(sort(otu_abd, decreasing = T))
  # Fetch hash otus names
  hash_otus<-taxa_names(physeq)
  # Associate hash to rank abundance prefix
  taxa_names(physeq)<-paste("OTU", match(hash_otus, rank_otus), sep = "_")
  # Create correspondence table for tracability
  correspondence<-data.frame(HashID = hash_otus,
                             RankID = taxa_names(physeq))
  # Write correspondence table to file
  write.csv(x = correspondence,file = hash2id_file,row.names = FALSE)
  return(physeq)
}
#' Convert taxa/ID through a correspondence table
#' @param physeq \code{\link{phyloseq-class}}
#' @param hashTable  A \code{\link{data.frame}} with at least a **UsearchID** and **HashID** columns.
#' @param ... The subsetting expression that should be applied to the
#'  correspondance table \code{hashTable}. This is passed on to \code{\link[base]{subset}}.
#'
#'
#'
#' @return physeq-renamed \code{\link{phyloseq-class}} with OTU renamed.
#' @export
#'
#' @examples
id2hash<-function(physeq, hashTable, ...){
  # Select only rows on interest, assembly method or min overlap for ex.
  # keep only ID columns
  intermediate.df<-subset(hashTable, subset = ..., select = c("UsearchID", "HashID"))
  # Drop unused levels
  intermediate.df<-droplevels(intermediate.df)
  # Replace rownames by USEARCH ID to speed up replacement.
  rownames(intermediate.df)<-intermediate.df$UsearchID
  # Discard unused column
  intermediate.df$UsearchID<-NULL

  # Fetch OTU ids
  otu_id<-taxa_names(physeq)
  # Replace USEARCH Id by Hash ID
  taxa_names(physeq)<-as.vector(intermediate.df[ otu_id, ])
  # Return modified phyloseq object
  return(physeq)
}
#' Taxa/OTU abundance in a given sample coupled with rank and taxonomic information
#'
#' @param physeq \code{\link{phyloseq-class}}
#' @param smp A single sample of interest. This integer or sample name will be
#' passed to \code{\link[phyloseq]{get_taxa}}.
#' @param level A character indicating the taxonomic rank of interest. Should be in
#' \code{\link[phyloseq]{rank_names}}. Default is "Species".
#'
#' @return An abundance sorted \code{\link{data.frame}} with the following components
#' \describe{
#' \item{\code{Sample}}{Sample of interest}
#' \item{\code{Abundance}}{Taxa/OTU abundance in the sample of interest}
#' \item{\code{Rank}}{Taxa/OTU rank in the dataset. The higher the rank, the less abundant the OTU.}
#' \item{\code{level}}{Taxonomic rank of interest. Note that the column is named after
#' the provided value from \code{level}.}
#' }
#' @export
#' @seealso \code{\link[phyloseq]{get_taxa}}
#' @examples
get_taxa_nomy<-function(physeq, smp, level="Species"){
  # Fetch OTU/taxa in samples
  taxa_list<-get_taxa(physeq, smp)
  # Subset list to taxa detected
  taxa_list<-taxa_list[ taxa_list > 0 ]
  df<-data.frame(
    Sample = smp,
    Abundance = taxa_list,
    Rank = match( names(taxa_list), names(sort(taxa_sums(physeq),decreasing = T))),
    TaxLevel = as(tax_table(physeq)[names(taxa_list),level],"matrix"))
  # Order by decreasing abundance
  df<-df[order(-df$Abundance),]
  # Include taxonomic level
  colnames(df)<-c('Sample','Abundance','Rank',level)
  return(df)
}

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

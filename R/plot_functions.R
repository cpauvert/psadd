
## Plot
#' Pie Plot Taxonomy Assignations from Two Abundance Vector
#'
#' Displays an overlay of two pie chart: (i) inner and (ii) outer.
#' It was primarily developed to illustrate taxonomic assignations at both
#' levels at the same time
#' (e.g. \href{https://en.wikipedia.org/wiki/Phylum}{Phylum} and
#'  \href{https://en.wikipedia.org/wiki/Class_(biology)}{Class})
#'
#' @param inner_abd A \code{\link{numeric-class}} describing the inner pie chart
#' @param outer_abd A \code{\link{numeric-class}} describing the outer pie chart
#' @param title A \code{\link{character}} specifying the plot title
#' @param COL.TAXO A palette (\code{\link{list}}) containing (at least) two named elements
#' @param EDG Number of edges
#'
#' @details \code{COL.TAXO} should contain at least two named elements \code{DIV} and \code{CLASS}.
#' These vectors describe the inner and outer pie chart color, respectively.
#'
#' @author Thomas Fort
#' @references See Fort et al. (2016) in \href{https://doi.org/10.7717/peerj.2656}{PeerJ}.
#'
#' @seealso \code{\link{floating.pie}} (in \code{\link{plotrix}} package)
#'
#' @importFrom plotrix floating.pie
#' @importFrom graphics par pie
#'
#' @return Nothing
#' @export
#'
plot_pie_vector<-function(inner_abd,outer_abd, title , COL.TAXO = COL.TAXO, EDG = 200){
  pie(1, radius=1, init.angle=90, col=c('white'), border = NA, labels='', main=title)
  plotrix::floating.pie(0,0,outer_abd, radius=1, startpos=pi/2, col=COL.TAXO$CLASS,border=NA, edges=EDG)
  plotrix::floating.pie(0,0, 1,radius=0.52, col=c('white'), border = NA, edges=EDG)
  plotrix::floating.pie(0,0, inner_abd, radius=0.5, startpos=pi/2, col=COL.TAXO$DIV, border = NA, edges=EDG)
}

#' Pie Plot Taxonomy Assignations from Phyloseq
#'
#' Displays an overlay of two pie chart: (i) inner (e.g. \href{https://en.wikipedia.org/wiki/Phylum}{Phylum})
#' and (ii) outer (e.g. \href{https://en.wikipedia.org/wiki/Class_(biology)}{Class})
#'
#'
#' @param physeq \code{\link{phyloseq-class}} with a \code{\link{sample_data-class}}
#' @param title. A \code{\link{character}} specifying the plot title
#' @param condition A \code{\link{character}} value of the "Description" column in \code{\link{sample_data-class}}
#' @param COL.TAXO.  COL.TAXO A palette (\code{\link{list}}) containing (at least) two named elements
#'
#' @details It was primarily developed to illustrate taxonomic assignations at both
#' levels at the same time.
#'
#' @return A \code{\link{list}} of two elements (\code{inner} and \code{outer}) containing taxa names.
#' @author Charlie PAUVERT
#' @references Adapted from Fort et al. (2016) in \href{https://doi.org/10.7717/peerj.2656}{PeerJ}.
#' @seealso \code{\link{plot_pie_vector}}, \code{\link{psmelt}}
#' @import phyloseq
#' @importFrom stats aggregate
#' @export
#'
plot_pie<-function(physeq, title.,condition,COL.TAXO. = COL.TAXO){
  # Extract sample data.frame and merge with taxonomy information
  # The resulting data.frame is decomposed to have a
  df<-psmelt(physeq)

  inner_tax<-aggregate(Abundance ~ Phylum+Description, data=df,sum,simplify=T)
  outer_tax<-aggregate(Abundance ~ Class+Description, data=df,sum,simplify=T)

  plot_pie_vector(subset(inner_tax,Description==condition, select = "Abundance",drop= TRUE),
                  subset(outer_tax,Description==condition, select = "Abundance",drop= TRUE),
                  title = title.,
                  COL.TAXO = COL.TAXO., EDG = 200)
  list_taxa<-list(
    inner=subset(inner_tax,Description==condition, select = 1,drop= TRUE),
    outer=subset(outer_tax,Description==condition, select = 1,drop= TRUE))
  return(list_taxa)
}

#' Plot Taxonomy Assignations from Phyloseq with legend.
#'
#' @inheritParams plot_pie
#' @seealso \code{\link{plot_pie}}
#'
#' @import phyloseq
#' @importFrom graphics legend
#' @export
plot_pie_with_legend <- function(physeq, title., condition) {
  par(mfrow = c(1, 2))
  list_taxa <- plot_pie(physeq, title., condition)
  legend(x = "bottomleft", legend = list_taxa$inner, col = COL.TAXO$DIV,
         bty = "n", pch = 15, ncol = 1, border = 0, pt.cex = 2, cex = 0.8,
         title = "Phylum")
  pie(1, radius = 1, init.angle = 90, col = c("white"), border = NA,
      labels = "", main = "Taxonomic Assignations Legend")
  legend(x = "top", legend = list_taxa$outer, col = COL.TAXO$CLASS, bty = "n",
         pch = 15, ncol = 1, border = 0, pt.cex = 2, cex = 0.8, title = "Class")
}
#' Plot Mean Variance OTU in order to evaluate overdispersion
#'
#' The \code{\link{mean}} and variance (\code{\link{var}}) of \code{physeq} OTU table are computed.
#' Overdispersion can be visually assessed by the generated \code{\link{plot}}. The putative linear relation
#' between mean and variance is added with \code{\link{abline}} to guide the analysis.
#'
#' @param physeq \code{\link{phyloseq-class}}
#' @param log A \code{\link{logical}} stating whether to apply \code{\link{log1p}} transform to both axis.
#' @param title A title to the mean-variance plot
#'
#' @details The \code{\link{log1p}} transform was chosen to limit \code{\link{NaN}} values with null values.
#'
#' @references Modified from \href{https://doi.org/10.1186/s40168-016-0208-8}{Thorsen and Brejnrod et al. (2016)}.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object
#' @seealso \code{\link{taxa_stats}}
#'
#' @export
#' @import ggplot2
#' @import phyloseq
#' @importFrom stats var
#' @examples
#' require(phyloseq)
#' data(GlobalPatterns)
#' (mv1<-plot_mv(GlobalPatterns, title = "Mean-Variance Plot with log-transform"))
#' mv1.nolog<-plot_mv(GlobalPatterns,log = FALSE, title = "Mean-Variance Plot without log-transform")
plot_mv<-function(physeq,log=TRUE,title = NULL){
  p<-ggplot(
    data.frame(
      m=taxa_stats(physeq,mean),
      v=taxa_stats(physeq,var)),
    aes(x = m,y = v)) +
    geom_point()+
    geom_abline(intercept = 0,slope = 1)+
    labs(x = "Mean", y = "Variance")
  if(log){
    p<-p+scale_y_continuous(trans = "log1p")+scale_x_continuous(trans = "log1p")
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}


#' Plot taxa/OTU table to evaluate sparsity.
#'
#'
#' Taxa/OTU are ordered by their prevalence in the dataset, and samples are
#' ordered by their richness (or completedness).
#'
#' @inheritParams plot_mv
#'
#' @references Adapted from \href{https://doi.org/10.1186/s40168-016-0208-8}{Thorsen and Brejnrod et al. (2016)}.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object
#' @seealso \code{\link{plot_mv}},\code{\link{taxa_prev}}, \code{\link{sample_prev}}
#' @export
#' @import phyloseq
#' @import ggplot2
#' @examples
#' \dontrun{
#' require(phyloseq)
#' data(enterotype)
#' (entero.sparse<-plot_sparsity(enterotype, title = "Sparsity Plot")
#' }
plot_sparsity<-function(physeq, title = NULL){
  # Melt the OTU table and merge associated metadata
  df<-psmelt(physeq)[,1:3]# keep OTU, Sample, Abundance columns
  # Order OTU by their occurrence/prevalence
  df$OTU<-factor(x = df$OTU,
                 levels = names(sort(taxa_prev(physeq))) )
  # Order Samples by their richness (or completedness)
  df$Sample<-factor(x = df$Sample,
                    levels = names(sort(sample_prev(physeq))) )
  # Transform Abundance numeric vector into boolean presence/absence
  df$Abundance<-df$Abundance > 0
  # Plot OTU table
  p<-ggplot(df[df$Abundance,],
            aes(x=Sample, y=OTU))+
    geom_tile()+scale_x_discrete(drop = F) + scale_y_discrete(drop=F) +
    theme_minimal()+theme(axis.text.x = element_blank(),
                          axis.text.y = element_blank(),
                          axis.ticks = element_blank())+
    labs(x= paste(length(levels(df$Sample)),"samples - from low to high richness"),
         y= paste(length(levels(df$OTU)),"OTU -- from rare to ubiquist"))
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}

#' Interactive Taxonomy plot with Krona from a phyloseq object
#'
#' Construct and run a Krona Chart to compare taxonomic assignations between
#' different conditions.
#'
#' @param physeq \code{\link{phyloseq-class}} with a \code{\link{taxonomyTable-class}}
#' @param output A \code{\link{character}} stating the output filename for Krona Chart and the directory in which Krona files will be created.
#' @param variable A \code{\link{character}} indicating which variable from sample data to use.
#' @param trim Should spaces and brackets be converted to underscore automatically.
#'
#' @import phyloseq
#' @importFrom utils browseURL write.table
#' @export
#'
#' @references \href{https://github.com/marbl/Krona/wiki/KronaTools}{KronaTools}.
#' @examples
#' \dontrun{
#' require(phyloseq)
#' data(GlobalPatterns)
#' plot_krona(GlobalPatterns,"GP-krona", "SampleType")# issues with brackets
#' plot_krona(GlobalPatterns,"GP-krona", "SampleType",trim=T)
#' }
plot_krona<-function(physeq,output,variable, trim=F){
  # Check if KronaTools are installed.
  if( system(command = "which ktImportText",
              intern = FALSE,
              ignore.stdout = TRUE)) {
  stop("KronaTools are not installed. Please see https://github.com/marbl/Krona/wiki/KronaTools.")
  }
  if( is.null(tax_table(physeq)) ){
    stop("No taxonomy table available.")
  }
  if( ! variable %in% colnames(sample_data(physeq))){
    stop(paste(variable, "is not a variable in the sample data."))
  }
  if (trim == FALSE) {
    spec.char<- grepl(" |\\(|\\)", as(sample_data(physeq),"data.frame")[,variable] )
    if(sum(spec.char > 0 )){
      message("The following lines contains spaces or brackets.")
      print(paste(which(spec.char)))
      stop("Use trim=TRUE to convert them automatically or convert manually before re-run")
    }
  }
  # Melt the OTU table and merge associated metadata
  df<-psmelt(physeq)
  # Fetch only Abundance, Description and taxonomic rank names columns
  df<-df[ ,c("Abundance", variable, rank_names(physeq)) ]
  # Make sure there are no spaces left
  df[,2]<-gsub(" |\\(|\\)","",df[,2])
  # Convert the field of interest as factor.
  df[,2]<-as.factor(df[,2])

  # Create a directory for krona files
  dir.create(output)

  # For each level of the Description variable
  # Abundance and taxonomic assignations for each OTU are fetched
  # and written to a file that would be processed by Krona.
  for( lvl in levels(df[,2])){
    write.table(
      df[which(df[, 2] == lvl & df[,1] != 0), -2]
      file = paste0(output,"/",lvl, "taxonomy.txt"),
      sep = "\t",row.names = F,col.names = F,na = "",quote = F)
  }
  # Arguments for Krona command
  # taxonomic file and their associated labels.
  krona_args<-paste(output,"/",levels(df[,2]),
                    "taxonomy.txt,",
                    levels(df[,2]),
                    sep = "", collapse = " ")
  # Add html suffix to output
  output<-paste(output,".html",sep = "")
  # Execute Krona command
  system(paste("ktImportText",
               krona_args,
               "-o", output,
               sep = " "))
  # Run the browser to visualise the output.
  browseURL(output)
}

#' Violin distribution plot
#'
#' Violin plot to visualise sample sequence number distribution.
#'
#' @inheritParams plot_sparsity
#' @param x A \code{\link{character}} describing the variable from
#' the \code{\link[phyloseq]{sample_data}} object to plot against.
#'
#' @import ggplot2
#' @import phyloseq
#' @export
#' @return A \code{\link[ggplot2]{ggplot}} object.
#' @examples
#' require(phyloseq)
#' data("soilrep")
#' head( sample_data(soilrep) )
#' plot_violin(physeq = soilrep, x = "warmed")
plot_violin<-function(physeq, x, title = NULL){
  # Compute sample sums
  smp_sums<-sample_sums(physeq)
  # Extract sample data if any.
  smp_df<-as( sample_data(physeq), "data.frame")
  # Merge both information: sample data and sequence number
  df<-merge(smp_df, data.frame(SequenceNumbers = smp_sums), by = "row.names")
  # Violin Plot horizontal
  p<-ggplot(df, aes_string(x,y="SequenceNumbers", fill=x))+
    geom_violin()+
    geom_jitter(alpha=0.5)+
    coord_flip()+
    theme_minimal()
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}

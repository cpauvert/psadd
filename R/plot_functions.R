
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

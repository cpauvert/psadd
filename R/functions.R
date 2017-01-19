
# Functions file


## Additions to phyloseq package
## Statistics
#' Apply statistics function to taxa abundance vector
#'
#' It generalize the \code{\link{taxa_sums}} function and enable the use of several summary functions.
#' @param physeq \code{\link{otu_table-class}}, or \code{\link{phyloseq-class}}
#' @param ... The function to be applied to the OTU table
#' @returns A \code{\link{numeric-class}} with length equal to the number of species
#'  in the table, name indicated the taxa ID, and value equal to the result of
#'  the function applied through \code{...}.
#' @example
#' data(esophagus)
#' taxa_stats(esophagus, mean)[1:6]
#' taxa_stats(esophagus, var)[1:6]
taxa_stats<-function (physeq,...) 
{
  m <- otu_table(physeq)
  marg<-ifelse(taxa_are_rows(m), yes = 1, no = 2)
  return(apply(m, marg, ...))
}

# Based on `taxa_sums`
# Return a named vector of taxa/OTU with the number of samples where
# its abundance is not 0
 
taxa_occ<-function (physeq) 
{
  m <- otu_table(physeq)
    if (taxa_are_rows(m)) {
      rowSums(m > 0)
    }
    else {
      colSums(m > 0)
    }
}

# Return a named vector of samples with the number of taxa/OTU 
# seen in each sample
sample_occ<-function (physeq) 
{
  m <- otu_table(physeq)
  if (taxa_are_rows(m)) {
    colSums(m > 0)
  }
  else {
    rowSums(m > 0)
  }
}

# Subset samples and prune OTU with zeros abundance 
subset_samples_no_zero<-function(physeq,...){
  x<-subset_samples(physeq,...)
  x<-prune_taxa( 
    taxa_sums(x)!=0 , # TRUE/FALSE OTU Vector to prune
    x )
  return(x)
}

# sample_specific_otu<-function(physeq,...){
#   # Fetch only the sample of interest. Ex: Ctrl.ex
#   x<-subset_samples_no_zero(physeq,...)
#   print(dim(otu_table(x)))
#   # original OTU table
#   m<-otu_table(physeq)
#   print(dim(m))
#   if(taxa_are_rows(physeq)){
#     # Subset only OTU in the specific sample
#     # AND on every other samples of the OTU Table
#     m<-m[ taxa_names(x), ! sample_names(physeq) %in% sample_names(x) ]
#   } else {
#     # Subset only OTU in the specific sample
#     # AND on every other samples of the OTU Table
#     m<-m[ ! sample_names(physeq) %in% sample_names(x), taxa_names(x) ]
#   }
#   return(m)
# }
# taxa_occ(rsb)[taxa_names(rsb_ctrl_ex)] - 1
relative_abundance<-function(physeq){
  # return relative abundance of OTU in phyloseq object
  abd<-taxa_sums(physeq)
  return(abd/sum(abd))
}

relative_abd_without_specific_sample<-function(physeq,specific){
  # Discard sample of interest from the global otu table
  m<-prune_samples(! sample_names(physeq) %in% sample_names(specific),physeq)
  otu_to_keep<-taxa_names(specific)
  rel_abd<-relative_abundance(m)
  return(rel_abd[ otu_to_keep ])
}


## Data frames
ctrl_samples_stats<-function(physeq,physeq_ctrl){
  # RA stands for relative abundance
  # RA of control otu
  RA_ctrl_otu<-relative_abundance(physeq_ctrl)
  RA_without_ctrl_sample<-relative_abd_without_specific_sample(physeq, physeq_ctrl)
  # Prevalence of control OTU in the global dataset, corrected by the number of control sample (1 usually)
  otu_prevalence<-taxa_occ(physeq)[taxa_names(physeq_ctrl)] - nsamples(physeq_ctrl)
  df<-data.frame(
    RA_In_Ctrl=RA_ctrl_otu,
    RA_In_All=RA_without_ctrl_sample,
    Ratio=RA_ctrl_otu/RA_without_ctrl_sample,
    Prevalence=otu_prevalence)
  return(df)
}

## Plot
plot_pie_vector<-function(inner_abd,outer_abd, title , COL.TAXO = COL.TAXO, EDG = 200){
  require(plotrix)
  pie(1, radius=1, init.angle=90, col=c('white'), border = NA, labels='', main=title)
  floating.pie(0,0,outer_abd, radius=1, startpos=pi/2, col=COL.TAXO$CLASS,border=NA, edges=EDG)
  floating.pie(0,0, 1,radius=0.52, col=c('white'), border = NA, edges=EDG)
  floating.pie(0,0, inner_abd, radius=0.5, startpos=pi/2, col=COL.TAXO$DIV, border = NA, edges=EDG)
}

plot_pie<-function(physeq, title.,condition,COL.TAXO. = COL.TAXO){
  require(phyloseq)
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

# plot_pie(rsb_wo_ctrl_otu,"Organic--Raw","Bio")
# plot_pie(rsb_wo_ctrl_otu_tr,"Organic -- Transform", "Bio")

# Plot Mean Variance OTU in order to evaluate overdispersion
plot_mv<-function(physeq,...){
  plot( 
    x = taxa_stats(physeq, mean),
    y = taxa_stats(physeq, var),
    ...)
  abline(a = 1,b = 0)
}

# Plot sparsity matrix
## based on Thorsen and Brejnrod et al. 2016
plot_sparsity<-function(physeq){
  # Melt the OTU table and merge associated metadata
  df<-psmelt(physeq)[,1:3]# keep OTU, Sample, Abundance columns
  # Order OTU by their occurrence/prevalence
  df$OTU<-factor(x = df$OTU,
                 levels = names(sort(taxa_occ(physeq))) )
  # Order Samples by their richness (or completedness)
  df$Sample<-factor(x = df$Sample,
                    levels = names(sort(sample_occ(physeq))) )
  p<-ggplot(df[df$Abundance,],
         aes(x=Sample, y=OTU))+
    geom_tile()+scale_x_discrete(drop = F) + scale_y_discrete(drop=F) +
    theme_minimal()+theme(axis.text.x = element_blank(),
                          axis.text.y = element_blank(),
                          axis.ticks = element_blank())+
    labs(x= paste(length(levels(df$Sample)),"samples - from low to high richness"),
         y= paste(length(levels(df$OTU)),"OTU -- from rare to ubiquist"),
         title="Sparsity Plot - based on Thorsen and Brejnrod et al. 2016")
  return(p)
}

# Miscellaneous

COL.TAXO<-list(
  DIV = c("#2b61ab", "#ed3247", "#259831", "#e87411", 
          "#199e81", "grey"),
  ASCO = c("#cd5b5b", "#00ced2", "#ff6525", 
           "#000083", "#aefe2d", "#a0522c", "#dc123a", "#01415d", "#007900", 
           "#bb58cf", "grey"),
  BASIDIO = c("#c61484", "#17176f", "#f3e38b", 
              "#f1817f", "#2e4e4d", "#ff4700", "#9acd30", "#2191ff", "#9e0000", 
              "#6b59cf", "#b2f0ef", "grey"),
  CLASS = c("#cd5b5b", "#00ced2", 
            "#ff6525", "#000083", "#aefe2d", "#a0522c", "#dc123a", "#01415d", 
            "#007900", "#bb58cf", "grey", "#c61484", "#17176f", "#f3e38b", 
            "#f1817f", "#2e4e4d", "#ff4700", "#9acd30", "#2191ff", "#9e0000", 
            "#6b59cf", "#b2f0ef", "grey", "grey", "grey", "grey", "grey"))

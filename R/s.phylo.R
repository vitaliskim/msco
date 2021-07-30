#' Species phylogeny generator
#'
#' Given species-by-site matrix (community), this function:
#' * uses the \link[taxize]{tax_name} function to obtain (from the
#'  [NCBI](https://www.ncbi.nlm.nih.gov/) or [ITIS](https://www.itis.gov/) online
#'   databases) the genus and family taxa levels of species in the community.
#'    If NCBI is used, getting an API key is recommended. See \link[taxize]{tax_name}
#'     for more information. NCBI is used as default in this function;
#' * uses the `phylo.maker` function to obtain the phylogeny (an object
#'  of `class`: "\code{phylo}") of species in the community using taxa obtained above;
#' * computes the phylogenetic distance matrix using the \link[ape]{cophenetic.phylo}
#'  function and the phylogeny obtained above as input;
#' * plots the phylogenetic tree using the phylogeny obtained above.
#'
#' @param s.data A species-by-site presence/absence `data.frame` with entries indicating
#' occurrence (1) and non-occurrence (0) of species in a site. The rows should have species'
#'  scientific names following [binomial nomenclature](https://en.wikipedia.org/wiki/Binomial_nomenclature),
#'   with no initials.
#' @param database The online database used to obtain the taxonomic names (species,
#'  genus and family) for a given rank (species list in this function). The options
#'   are "ncbi" or "itis". We used "ncbi" as default in this function.
#' @param taxa.levels Species taxa (i.e. a `data.frame` with species, genus, and family as `colnames`) used
#'  in extracting phylogenetic distance matrix between species. If supplied, `taxa.levels` won't be
#'   computed from online repositories. Taxa provision is highly recommended.
#' @param Obs.data A logical indicating if `s.data` should be output.
#' @param obs.taxa A logical indicating if `taxa.levels` should be output.
#' @param phy.d.mat A logical indicating if phylogenetic distance matrix should be output.
#' @param phylo.plot Boolean value indicating if phylogenetic tree (cluster dendrogram) should be plotted.
#' @return Returns a `list` with the following outputs:
#' * `s.data`: &nbsp;A `data.frame` with sites as columns and species as rows.
#' * `taxa.levels`: &nbsp;A `data.frame` with the following columns:
#'     + `species`: &nbsp;Species names in `s.data`.
#'     + `genus`: &nbsp;Genus names of species in `s.data`.
#'     + `family`: &nbsp;Family names of species in `s.data`.
#' * `p.d.matrix`: &nbsp;A symmetric `matrix` with `dimnames` as species and entries indicating the
#'  phylogenetic distance between any two of them (species).
#' * `phylo.plot`: &nbsp;A phylogenetic tree (cluster dendrogram) of species in `s.data`
#'
#' @references
#' \enumerate{
#'  \item{Binomial nomenclature' (2020) *Wikipedia*. Available at:
#'  <https://en.wikipedia.org/wiki/Binomial_nomenclature> (Accessed: 09 November 2020).}
#'
#'  \item{Lagat, V. K., Latombe, G. and Hui, C., 202X. *Dissecting the effects of
#'   neutral encounter versus functional traits on multi-order species interactions
#'    and co-occurrence with generalised B-spline modelling*. Upcoming.}
#' }
#' @examples
#' \dontrun{
#'
#' remotes::install_github("jinyizju/V.PhyloMaker", force = TRUE)
#' library(V.PhyloMaker)
#' my.path <- system.file("extdata/gsmdat", package = "msco")
#' setwd(my.path)
#' s.data <- get(load("s.data.csv"))
#' taxa <- get(load("taxa.levels.csv"))
#'
#' my.s.phylo <- msco::s.phylo(s.data, database = "ncbi", obs.taxa=TRUE,
#'  taxa.levels = taxa, Obs.data=TRUE, phy.d.mat=TRUE, phylo.plot = TRUE)
#'
#' my.s.data <- my.s.phylo$s.data
#' my.s.data
#'
#' my.taxa <- my.s.phylo$taxa.levels
#' my.taxa
#'
#' my.p.d.mat <- my.s.phylo$phylogenetic.distance.matrix
#' my.p.d.mat
#'
#' }
#' @export
#' @md
s.phylo <- function(s.data, database = "ncbi", obs.taxa=TRUE, taxa.levels = NULL, Obs.data=TRUE,
                   phy.d.mat=TRUE,  phylo.plot = TRUE){

  #taxa
  if(is.null(taxa.levels)){
    vee <- taxize::tax_name(sci = c(row.names(s.data)), get = c("genus","family"), db = "ncbi")
    if(database == "itis"){
      vee <- taxize::tax_name(sci = c(row.names(s.data)), get = c("genus","family"), db = "itis")
    }
    taxa <- vee[,c(2,3,4)]
    names(taxa)[1] <- "species"
    taxa <- taxa[stats::complete.cases(taxa),]
  }else {
    taxa <- taxa.levels
  }

  ##Phylogenetic distance matrix
  kimm <- import::here("V.PhyloMaker", "nodes.info.1", "GBOTB.extended", "phylo.maker")
  vdat <- kimm$phylo.maker(sp.list = taxa, tree = kimm$GBOTB.extended)$scenario.3
  p.d.mat <- as.matrix(ape::cophenetic.phylo(vdat))

  ##Create a list of all outputs
  phylo.vee <- list()
  if(nrow(p.d.mat) != nrow(s.data)){
    warning(paste(deparse(setdiff(base::chartr(" ", "_", row.names(s.data)), row.names(p.d.mat))),
                  " species failed to be binded to the tree", ". It was removed from s.data.",
                  sep = ""))
    st.data <- `row.names<-`(s.data, base::chartr(" ", "_", row.names(s.data)))
    stt.data <- st.data[intersect(row.names(st.data), row.names(p.d.mat)),]
    sttt.data <- `row.names<-`(stt.data, base::sub("_", " ", row.names(stt.data)))
    if(Obs.data==TRUE){
      phylo.vee$s.data <- as.data.frame(sttt.data)
    }
  }else if(nrow(p.d.mat) == nrow(s.data)){
    st.data <- `row.names<-`(s.data, base::chartr(" ", "_", row.names(s.data)))
    stt.data <- st.data[intersect(row.names(st.data), row.names(p.d.mat)),]
    sttt.data <- `row.names<-`(stt.data, base::sub("_", " ", row.names(stt.data)))
    if(Obs.data==TRUE){
      phylo.vee$s.data <- as.data.frame(sttt.data)
    }
  }
  if(obs.taxa==TRUE){
    phylo.vee$taxa.levels <- taxa
  }

  if(phy.d.mat==TRUE){
    phylo.vee$phylogenetic.distance.matrix <- as.data.frame(p.d.mat)
  }


  if(phylo.plot == TRUE){
    grDevices::pdf(file = paste0(system.file("ms", package = "msco"), "/Phylogenetic.tree.pdf"), paper="a4r", height = 8.27, width = 11.69)
    hc <- stats::as.hclust(vdat)
    dend <- stats::as.dendrogram(hc)
    graphics::par(mar = c(12, 6, 4, 0)) # leave space for the labels
    graphics::plot(dend, ylab = "Phylogenetic distance", main = "Cluster Dendrogram")
    # phylo.vee$phylo.plot <- grDevices::recordPlot()
    grDevices::dev.off()
    base::system(paste0('open "', paste0(system.file("ms", package = "msco"), "/Phylogenetic.tree.pdf"), '"'))
  }

  return(phylo.vee)
}



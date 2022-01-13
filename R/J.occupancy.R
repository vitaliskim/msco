
#' A schematic figure of the archetypes
#'
#' A schematic diagram illustrating nine possible archetypes (from the null model test) of the
#'  patterns of species co-occurrences in ecological communities. The archetypes are denoted
#'  \eqn{\{Ai: i \in (1:9)\}}. See the details below.
#'
#' @return The `Arch_schem` function returns a **schematic diagram** of the archetypes of species
#'  co-occurrence patterns (denoted by \eqn{\{Ai: i \in (1:9)\}}), with the following components:
#'  \item{**Archetype**}{**Description/Interpretation**}
#'  \item{A1}{The joint occupancy value of the observed community matrix (observed; dark solid line) is above the null
#'   model. This means the null hypothesis (i.e. a statement that imply any change in the observed patterns do
#'    not reflect any community assembly process as underlying cause) should be rejected, confirming the presence of a
#'    mechanism of interest being tested (Lagat *et al.,* 2021a). It is typical of a community whose
#'    species are positively associated (or aggregated) more often than would be expected by chance.
#'     Such patterns of community structure may arise from a number of ecological mechanisms
#'      including environmental filtering or shared habitat requirements (Cordero and Jackson,
#'       2019).}
#'  \item{A2}{The observed is greater than null expectation for `i = 2` but within null expectation
#'   for \eqn{i \ge 3}. This implies a pairwise metric detects a non-random pattern of the community
#'   structure, but when higher order species are considered, a random pattern is produced. This is
#'    typical of a community whose species are aggregated more often than by chance in sites with few
#'     species than in sites with many species (Lagat *et al.,* 2021a).}
#'  \item{A3}{The observed is greater than null expectation for lower orders, within null expectation for
#'   medium orders, and less than null expectation for higher orders. This means species co-occur more
#'    often than by chance in sites with few species, but are segregated more often than by chance in
#'     sites with many species, depicting a community structured by two different community assembly
#'      processes (Lagat *et al.,* 2021a).}
#'  \item{A4}{The observed is within null expectation for `i = 2` but greater than null expectation for
#'   \eqn{i \ge 3}. This means when use pairwise co-occurrence is used, the null hypothesis is not
#'    rejected, but when joint occupancy is used, the same null hypothesis is rejected. I.e., pairwise
#'     co-occurrence fails at detecting patterns of aggregation for sites with many species, i.e. a
#'      type II error or false negative (Lagat *et al.,* 2021a).}
#'  \item{A5}{The observed is within the null expectation for all orders \eqn{i \ge 2}, implying the
#'   test is not statistically significant. This has been ecologically inferred to mean ecological
#'    communities are random and that no community assembly processes or mechanisms influence their
#'     structure (Lagat *et al.,* 2021a; Cordero and Jackson, 2019; Gotelli and Sounding, 2001).}
#'  \item{A6}{The observed is within the null expectation for `i = 2` but less than the null expectation
#'   for \eqn{i \ge 3}. This means when use pairwise co-occurrence is used, the null hypothesis is not
#'    rejected, but when joint occupancy is used, the same null hypothesis is rejected. I.e., pairwise
#'     co-occurrence fails at detecting patterns of segregation for sites with many species, i.e. a
#'      type II error or false negative (Lagat *et al.,* 2021a).}
#'  \item{A7}{The observed is less than null expectation for lower orders, within null expectation
#'   for medium orders, and greater than null expectation for higher orders. Implying species are
#'    segregated more often than would be expected by chance in sites with few species, but co-occur
#'     more often than by chance in sites with many species, depicting a community structured by
#'      two different community assembly processes (Lagat *et al.,* 2021a).}
#'  \item{A8}{The observed is less than null expectation for `i = 2` but within null expectation for
#'   \eqn{i \ge 3}. This means a pairwise metric detects a non-random pattern of the community structure,
#'   but when higher order species are considered, a random pattern is produced. This is typical of a
#'    community whose species are segregated more often than by chance in sites with few species than
#'     in sites with many species (Lagat *et al.,* 2021a).}
#'  \item{A9}{The joint occupancy value of the observed community matrix (dark solid line) is below the
#'   null model. This means the null hypothesis should be rejected, confirming the presence of a
#'    mechanism of interest being tested (Lagat *et al.,* 2021a). It is typical of a community
#'     structured by inter-specific competition or limiting similarity, though predation might also
#'      generate similar patterns (Hein et al. 2014).}
#' @references
#' \enumerate{
#'
#'  \item{Cordero, R.D. and Jackson, D.A. (2019). Species-pair associations, null models, and tests of
#'    mechanisms structuring ecological communities. *Ecosphere* **10.** <https://doi.org/10.1002/ecs2.2797>}
#'
#'  \item{Gotelli, N. J. and Sounding, E. (2001). Research frontiers in null model analysis. *Glob. Ecol.
#'    Biogeogr.* **10**, 337-343. <https://doi.org/10.1046/j.1466-822X.2001.00249.x>}
#'
#'  \item{Hein et al. (2014). Fish introductions reveal the temperature
#'   dependence of species interactions. *Proc. R. Soc. B Biol. Sci.* **281**.
#'    <https://doi.org/10.1098/rspb.2013.2641>}
#'
#'  \item{Lagat, V. K., Latombe, G. and Hui, C. (2021a). *A multi-species co-occurrence index to
#'   avoid type II errors in null model testing*. DOI: `<To be added>`.}
#' }
#' @note `Arch_schem` is not a generic function which can take in any dataset and give the outputs,
#'  but a path to a schematic diagram saved in this package. A representational figure from empirical,
#'   simulated or any known `.csv` binary data matrices can be accessed with \link[msco]{Jo.plots} function.
#'
#' @export
#' @md

Arch_schem <- function(){
  print(noquote("Check msco's 'schem' folder in your R version's directory for a 'Schematic_archetype_figures.pdf' file."))
}

#' Joint occupancy parametric and null model plots
#'
#' Plots the null model and joint occupancy decline with order (number of species)
#' and fits the decline to exponential, power law and exponential-power law
#' parametric models, respectively.
#'
#' @param jo_Obj A joint occupancy model object returned by the function \link[msco]{Jo.eng}.
#' @return Produces a figure consisting of the following plots:
#'  \item{(a)}{Joint occupancy decline.}
#'  \item{(b)}{Exponential regression of the joint occupancy decline.}
#'  \item{(c)}{Power law regression of the joint occupancy decline.}
#'  \item{(d)}{Exponential-power law regression of the joint occupancy decline.}
#'  \item{(e)}{Null model test.}
#' @details This function provides a visualization of the forms of joint
#'   occupancy decline and null model test. It offers information on:
#'  \itemize{
#'  \item{the outcomes of the null model test (through the appended
#'    archetype value on (e) plot) and}
#'  \item{the comparisons between the forms of joint occupancy decline (through
#'      the affixed `AIC` and `rsq` values on (b), (c) and (d) plots,
#'       respectively)}.
#'   }
#' @references Lagat, V. K., Latombe, G. and Hui, C. (2021a). *A multi-species co-occurrence
#'  index to avoid type II errors in null model testing*. DOI: `<To be added>`.
#' @examples
#' \dontrun{
#'
#' ex.data <- read.csv(system.file("extdata", "274.csv", package = "msco"))
#' jo_Obj <- msco::Jo.eng(ex.data, nReps = 999, All.plots = TRUE, s.dplot = FALSE, dig = 3)
#' jplots <- msco::Jo.plots(jo_Obj)
#' jplots
#'
#' ex.data2 <- read.csv(system.file("extdata", "22.csv", package = "msco"))
#' jo_Obj2 <- msco::Jo.eng(ex.data2, nReps = 999, All.plots = TRUE, s.dplot = FALSE, dig = 3)
#' jplots2 <- msco::Jo.plots(jo_Obj2)
#' jplots2
#'
#' ex.data3 <- read.csv(system.file("extdata", "78.csv", package = "msco"))
#' jo_Obj3 <- msco::Jo.eng(ex.data3, nReps = 999, All.plots = TRUE, s.dplot = FALSE, dig = 3)
#' jplots3 <- msco::Jo.plots(jo_Obj3)
#' jplots3
#'
#' ex.data4 <- read.csv(system.file("extdata", "65.csv", package = "msco"))
#' jo_Obj4 <- msco::Jo.eng(ex.data4, nReps = 999, All.plots = TRUE, s.dplot = FALSE, dig = 3)
#' jplots4 <- msco::Jo.plots(jo_Obj4)
#' jplots4
#'
#' }
#' @export
#' @md

Jo.plots <- function(jo_Obj){
  jo_Obj$all.plots
}

#' Expected value of joint occupancy for order \eqn{i} and its standard deviation
#'
#' This function computes joint occupancy (the average number of sites harbouring a given number of \eqn{i} species
#'  simultaneously), and its standard deviation.
#'
#' @param s.data A species-by-site presence/absence matrix with entries indicating
#' occurrence (1) and non-occurrence (0) of species in a site.
#' @param order Specific number of species for which joint occupancy and its standard deviation
#'  is computed.
#' @param metric The type of rescaling applied to the joint occupancy metric. Available options are:
#'  `Simpson_eqn` for Simpson equivalent, `Sorensen_eqn` for Sorensen equivalent, and `raw` for the
#'   raw form of index without rescaling.
#' @return Returns a `list` with the following outputs:
#' \item{jo.val}{Joint occupancy value.}
#' \item{jo.sd}{The standard deviation of `jo.val`.}
#' @references Lagat, V. K., Latombe, G. and Hui, C. (2021a). *A multi-species co-occurrence
#'  index to avoid type II errors in null model testing*. DOI: `<To be added>`.
#' @examples
#' ex.data <- read.csv(system.file("extdata", "274.csv", package = "msco"))
#' jo <- msco::j.occ(ex.data, order = 3, metric = "raw")
#' jo
#'
#' ex.data2 <- read.csv(system.file("extdata", "65.csv", package = "msco"))
#' jo2 <- msco::j.occ(ex.data2, order = 3, metric = "raw")
#' jo2
#' @export
#' @md

j.occ<-function(s.data, order, metric = "raw"){

  s.data <- as.matrix(s.data)
  s.data <- s.data[rowSums(s.data) > 0, ] ## Remove rows with no species
  richness <- colSums(s.data)
  p <- exp(lchoose(richness, order) - lchoose(nrow(s.data),order))
  jo.val <- sum(p)

  similarity_mat <- t(s.data) %*% s.data
  covmat <- exp(lchoose(similarity_mat, order) - lchoose(nrow(s.data), order))
  for (j in 1:ncol(s.data)) {
    for (k in 1:ncol(s.data)) {
      covmat[j, k] <- covmat[j, k] - p[j] * p[k]
    }
  }
  jo.var <- choose(nrow(s.data), order)/
    (choose(nrow(s.data), order) - 1) * sum(covmat)

  jo.sd <- sqrt(jo.var)

  jo.index <- list()
  if(metric == "raw"){
    jo.index$jo.val <- jo.val
  }else if(metric=="Simpson_eqn"){
    jo.index$jo.val <- jo.val/min(rowSums(s.data))
  }else if(metric=="Sorensen_eqn"){
    jo.index$jo.val <- jo.val/mean(rowSums(s.data))
  }else if((metric %in% c("raw", "Simpson_eqn", "Sorensen_eqn"))!=TRUE){
    stop("Wrong option for the joint occupancy metric provided. It must either be 'raw', 'Simpson_eqn', or 'Sorensen_eqn'.")
  }
  jo.index$jo.sd <- jo.sd

  return(jo.index)

}

#' Expected value of joint occupancy and its standard deviation for a range of orders
#'
#' This function computes joint occupancy (the average number of sites harbouring a
#'  given number of \eqn{i} species simultaneously) and its standard deviation for a
#'   range of orders (number of species).
#'
#' @param s.data A species-by-site presence/absence matrix with entries indicating
#' occurrence (1) and non-occurrence (0) of species in a site.
#' @param orders Range number of species for which joint occupancy and its standard
#'  deviation is computed.
#' @param metric The type of rescaling applied to the joint occupancy metric. Available options are:
#'  `Simpson_eqn` for Simpson equivalent, `Sorensen_eqn` for Sorensen equivalent, and `raw` for the
#'   raw form of index without rescaling.
#' @return Returns a `list` with the following outputs:
#' \item{jo.vals}{A vector of joint occupancy values for a range number of species (in `orders`).}
#' \item{jo.sds}{A vector of standard deviations of `jo.vals`.}
#' @references Lagat, V. K., Latombe, G. and Hui, C. (2021a). *A multi-species co-occurrence
#'  index to avoid type II errors in null model testing*. DOI: `<To be added>`.
#' @examples
#' ex.data <- read.csv(system.file("extdata", "274.csv", package = "msco"))
#' jos <- msco::j.occs(ex.data, orders = 1:nrow(ex.data), metric = "raw")
#' jos
#'
#' ex.data2 <- read.csv(system.file("extdata", "65.csv", package = "msco"))
#' jos2 <- msco::j.occs(ex.data2, orders = 1:nrow(ex.data), metric = "raw")
#' jos2
#' @export
#' @md

j.occs<-function(s.data, orders = 1:nrow(s.data), metric = "raw"){

  jo=0
  SDs=0
  for (i in orders) {
    jo[i]=j.occ(s.data, order = i, metric = metric)$jo.val
    SDs[i]=j.occ(s.data, order = i, metric = metric)$jo.sd
  }
  jo.inds <- list()
  jo.inds$jo.vals <- jo
  jo.inds$jo.sds <- SDs
  return(jo.inds)
}

#' Results on joint occupancy index (presented in Lagat et al., 2021a)
#'
#' This function allows the replication of the results presented in Lagat *et al*. (2021a). Executing
#'   `Jo.res()` therefore gives these outputs that are saved as `.RDS` files in `msco`. If the codes that
#'    produced these (saved) outcomes are desired, the codes below are made available.
#'
#' @note The function \link[msco]{Jo.res} is not for general use. We included it in this package to help
#'  the readers of Lagat *et al*. (2021a) paper, who may want to get a deeper understanding of how the results
#'   presented in this paper were arrived at. It also allows deeper scrutiny of Lagat *et al*. (2021a)'s
#'    findings.
#'
#' @return Returns all the results presented in Lagat *et al*. (2021a). To replicate these results,
#'  execute the following code:
#'
#'  ```
#'    RNGkind(sample.kind = "Rejection")
#'    set.seed(39)
#'    my.path <- system.file("extdata/myCSVs", package = "msco")
#'    setwd(my.path)
#'    my.files <- gtools::mixedsort(list.files(path = my.path, pattern = "*.csv"))
#'    Lag.res <- msco::mJo.eng(my.files,
#'                     algo = "sim2",
#'                     metric = "raw",
#'                     Archetypes = FALSE,
#'                     AICs = FALSE,
#'                     params = FALSE,
#'                     my.r2 = FALSE,
#'                     my.r2.s = TRUE,
#'                     best.mod2 = TRUE,
#'                     best.mod3 = TRUE,
#'                     params_c.i = TRUE
#'                    );Lag.res
#'
#'  ```
#'
#'  **Caveat:** The above code can take approximately 10 minutes to execute. It took 10.39014 minutes to run
#'   (and output results) on a 64 bit system with 8 GB RAM and 3.60 GHz CPU.
#'
#'  * __Fig. 3__ can be replicated using:
#'  ```
#'   grDevices::dev.new()
#'   msco:::nullmod_archs()
#'
#'  ```
#'  * __Fig. S1__ can be replicated using:
#'  ```
#'  my.path <- system.file("extdata/myCSVs", package = "msco")
#'  setwd(my.path)
#'  my.files <- gtools::mixedsort(list.files(path = my.path, pattern = "*.csv"))
#'  grDevices::dev.new()
#'  msco:::richness.variances(my.files)
#'
#'  ```
#'
#' @references Lagat, V. K., Latombe, G. and Hui, C. (2021a). *A multi-species co-occurrence
#'  index to avoid type II errors in null model testing*. DOI: `<To be added>`.
#'
#' @examples
#' \dontrun{
#'
#' ms.res <- msco::Jo.res()
#' ms.res$r2.s
#' ms.res$best.mod2
#' ms.res$best.mod3
#' ms.res$params_c.i
#' }
#' @export
#' @md

Jo.res <- function(){
  res <- readRDS(system.file("ms", "jo.res.RDS", package = "msco"))
  print(nullmod_archs())
  return(res)
}

#' Results on generalised B-spline modelling (presented in Lagat et al., 2021b)
#'
#' This function allows the replication of the results on generalised B-spline modelling, presented
#'  in Lagat *et al*. (2021b). Executing `gbsm.res()` therefore gives these outputs that are saved
#'   as `.RDS` files in `msco`. If the codes that produced these (saved) outcomes are desired, the
#'    codes below are made available.
#'
#' @note The function \link[msco]{gbsm.res} is not for general use. We included it in this package to help
#'  the readers of Lagat *et al*. (2021b) paper, who may want to get a deeper understanding of how the results
#'   presented in this paper were arrived at. It also allows deeper scrutiny of Lagat *et al*. (2021b)'s
#'    findings.
#'
#' @return Returns all the results presented in Lagat *et al*. (2021b). To replicate
#'
#' * __Figs. 1__, __3__, __4__, __5__, and __Tables 1__ and __S1__, execute the following code:
#'
#'
#'   ```
#'
#'    my.path <- system.file("extdata/gsmdat", package = "msco")
#'    setwd(my.path)
#'    s.data <- get(load("s.data.csv")) ##Species-by-site matrix
#'    t.data <- get(load("t.data.csv")) ##Species-by-trait matrix
#'    p.d.mat <- get(load("p.d.mat.csv")) ##Species-by-species phylogenetic distance matrix
#'    RNGkind(sample.kind = "Rejection")
#'    set.seed(1)
#'    gb.res <- msco::gbsm_m.orders(s.data,
#'                t.data,
#'                p.d.mat,
#'                metric = "Simpson_eqn",
#'                gbsm.model,
#'                orders = c(2:5, 8, 10, 15),
#'                d.f = 4,
#'                degree = 3,
#'                n = 1000,
#'                k = 5,
#'                p = 0.8,
#'                type = "k-fold",
#'                scat.plots = TRUE,
#'                response.curves = TRUE,
#'                j.occs.distrbn = TRUE,
#'                mp.plots = TRUE,
#'                start.range=c(-0.1,0)
#'              )
#'
#'    gb.res$contbn_table$`order 3`  ## Table 1
#'    gb.res$model.validation.table  ## Table S1
#'
#'
#'   ```
#'
#' * __Figs. S1, 2,__ and __S2__, execute the following codes:
#'
#'      + __Fig. S1__:
#'      ```
#'        remotes::install_github("jinyizju/V.PhyloMaker", force = TRUE)
#'        library(V.PhyloMaker)
#'        my.path <- system.file("extdata/gsmdat", package = "msco")
#'        setwd(my.path)
#'        s.data <- get(load("s.data.csv")) ##Species-by-site matrix
#'        taxa <- get(load("taxa.levels.csv")) ##Species taxa
#'        my.phylo.plot <- msco::s.phylo(s.data,
#'                                  database = "ncbi",
#'                                  obs.taxa = FALSE,
#'                                  taxa.levels = taxa,
#'                                  Obs.data = FALSE,
#'                                  phy.d.mat = FALSE,
#'                                  phylo.plot = TRUE)
#'
#'      ```
#'
#'      + __Fig. 2__:
#'      ```
#'        my.path <- system.file("extdata/gsmdat", package = "msco")
#'        setwd(my.path)
#'        s.data <- get(load("s.data.csv")) ##Species-by-site matrix
#'        t.data <- get(load("t.data.csv")) ##Species-by-trait matrix
#'        p.d.mat <- get(load("p.d.mat.csv")) ##Species-by-species phylogenetic distance matrix
#'        my.gbsm <- msco::gbsm(s.data,
#'                          t.data,
#'                          p.d.mat,
#'                          metric = "Simpson_eqn",
#'                          gbsm.model,
#'                          d.f = 4,
#'                          order.jo = 3,
#'                          degree = 3,
#'                          n = 1000,
#'                          b.plots = TRUE,
#'                          bsplines = "single",
#'                          scat.plot = FALSE,
#'                          response.curves = FALSE,
#'                          leg = 1,
#'                          start.range=c(-0.1,0)
#'                        )
#'
#'      ```
#'      + __Fig. S2__:
#'      ```
#'        my.path <- system.file("extdata/gsmdat", package = "msco")
#'        setwd(my.path)
#'        s.data <- get(load("s.data.csv")) ##Species-by-site matrix
#'        t.data <- get(load("t.data.csv")) ##Species-by-trait matrix
#'        p.d.mat <- get(load("p.d.mat.csv")) ##Species-by-species phylogenetic distance matrix
#'        RNGkind(sample.kind = "Rejection")
#'        set.seed(4)
#'        pe <- msco::pred.error.bands(s.data,
#'                            t.data,
#'                            p.d.mat,
#'                            metric = "Simpson_eqn",
#'                            gbsm.model,
#'                            d.f = 4,
#'                            simm = 10,
#'                            orders = c(2:5, 8, 10, 15),
#'                            degree = 3,
#'                            n = 1000,
#'                            start.range=c(-0.1,0)
#'                          )
#'      ```
#'
#'  **Caveat:** The above codes can collectively take approximately 7 minutes to execute
#'   (with prediction uncertainty plot taking 6 minutes alone). It took 7.3895 minutes
#'    to run (and output results) on a 64 bit system with 8 GB RAM and 3.60 GHz CPU.
#'
#' @references Lagat, V. K., Latombe, G. and Hui, C. (2021b). *Dissecting the effects
#'  of random encounter versus functional trait mismatching on multi-species
#'   co-occurrence and interference with generalised B-spline modelling*. DOI: `<To be added>`.
#'
#' @examples
#' \dontrun{
#'
#' gbs.res <- msco::gbsm.res()
#' gbs.res$contbn_table$`order 3`
#' gbs.res$model.validation.table
#'
#' }
#' @export
#' @md

gbsm.res <- function(){
  gres <- readRDS(system.file("ms", "gbsm.res.RDS", package = "msco"))
  # system(paste0('open "', paste0(system.file("ms", package = "msco"), "/gbsm.plots.pdf"), '"'))
  return(gres)
}

#' Results on `msco` illustration (presented in Lagat et al., 2021c)
#'
#' This function allows the replication of the results on `msco` R package illustration paper presented
#'  in Lagat *et al*. (2021c). Executing `msco.res()` therefore gives these outputs that are saved
#'   as `.RDS` files in `msco`. If the codes that produced these (saved) outcomes are desired, the
#'    codes below are made available.
#'
#' @note The function \link[msco]{msco.res} is not for general use. We included it in this package to help
#'  the readers of Lagat *et al*. (2021c) paper, who may want to get a deeper understanding of how the results
#'   presented in this paper were arrived at. It also allows deeper scrutiny of Lagat *et al*. (2021c)'s
#'    findings, and broader understanding of the main functionalities of `msco` R package.
#'
#' @return Returns all the results presented in Lagat *et al*. (2021c). To replicate
#'
#' * __Figs. 1__, __2__ and __Table 2__, execute the following code:
#'
#'
#'   ```
#'
#'    RNGkind(sample.kind = "Rejection")
#'    set.seed(14)
#'    ex.data <- read.csv(system.file("extdata", "251.csv", package = "msco"))
#'    j.en <- msco::Jo.eng(ex.data,
#'                algo = "sim2",
#'                metric = "raw",
#'                nReps = 999,
#'                dig = 3,
#'                s.dplot = FALSE,
#'                All.plots = TRUE,
#'                Jo.coeff = TRUE,
#'                my.AIC = TRUE,
#'                my.rsq = TRUE,
#'                Exp_Reg = TRUE,
#'                P.law_Reg = TRUE,
#'                Exp_p.l_Reg = TRUE,
#'                Obs.data = FALSE,
#'                Sim.data = FALSE,
#'                Jo_val.sim = FALSE,
#'                lab = FALSE,
#'                leg = FALSE,
#'                C.I_Jo_val.sim = FALSE,
#'                Jo_val.obs = TRUE,
#'                Metric = TRUE,
#'                Algorithm = TRUE,
#'                S.order = TRUE,
#'                nmod_stats = TRUE,
#'                Pt_Arch_Vals = TRUE,
#'                Atype = TRUE,
#'                p.n.plot = TRUE,
#'                trans = FALSE,
#'                m.n.plot = FALSE)
#'
#'    j.en$jo.coeff ## Table 1
#'    j.en$AIC; j.en$r2 ## Table 2
#'    j.en$nmod_stats ## Table 3
#'    grDevices::dev.new()
#'    j.en$all.plots
#'
#'   ```
#'
#' * __Fig. 4__, execute the following code:
#'
#'   ```
#'    RNGkind(sample.kind = "Rejection")
#'    set.seed(14)
#'    grDevices::dev.new()
#'    msco:::nullmod_archs2()
#'
#'   ```
#'
#' * __Fig. 5__, execute the following code:
#'
#'   ```
#'
#'    my.path <- system.file("extdata/gsmdat", package = "msco")
#'    setwd(my.path)
#'    s.data <- get(load("s.data.csv")) #Species-by-site matrix
#'    t.data <- get(load("t.data.csv")) #Species-by-trait matrix
#'    p.d.mat <- get(load("p.d.mat.csv")) #Species-by-species phylogenetic distance matrix
#'    RNGkind(sample.kind = "Rejection")
#'    set.seed(1)
#'    gb.res <- msco::gbsm_m.orders(s.data,
#'                t.data,
#'                p.d.mat,
#'                metric = "Simpson_eqn",
#'                gbsm.model,
#'                orders = c(3:5, 8, 10, 15, 20),
#'                d.f = 4,
#'                degree = 3,
#'                n = 1000,
#'                k = 5,
#'                p = 0.8,
#'                type = "k-fold",
#'                scat.plots = FALSE,
#'                response.curves = TRUE,
#'                j.occs.distrbn = FALSE,
#'                mp.plots = FALSE,
#'                start.range=c(-0.1,0)
#'              )
#'
#'
#'   ```
#' @references
#' \enumerate{
#' \item{Lagat, V. K., Latombe, G. and Hui, C. (2021a). *A multi-species co-occurrence
#'  index to avoid type II errors in null model testing*. DOI: `<To be added>`.}
#'
#'  \item{Lagat, V. K., Latombe, G. and Hui, C. (2021b). *Dissecting the effects of random
#'   encounter versus functional trait mismatching on multi-species co-occurrence and
#'    interference with generalised B-spline modelling*. DOI: `<To be added>`.}
#'
#'   \item{Lagat, V. K., Latombe, G. and Hui, C. (2021c). *`msco`: an R software package
#'   for null model testing of multi-species interactions and interference with
#'    covariates*. DOI: `<To be added>`.}
#'
#' }
#'
#' @examples
#' \dontrun{
#'
#' ms.res <- msco::msco.res()
#' ms.res$nmod_stats ## Table 2
#'
#' }
#' @export
#' @md

msco.res <- function(){
  mres <- readRDS(system.file("ms", "msco.res.RDS", package = "msco"))
  # system(paste0('open "', paste0(system.file("ms", package = "msco"), "/msco.illus.plots.pdf"), '"'))
  # saveRDS(j.en, file = paste0(system.file("ms", package="msco"), "/msco.res.RDS"))

  return(mres)
}


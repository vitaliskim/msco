#' Joint occupancy model engine
#'
#' This function is the engine behind the null model testing of species co-occurrence patterns,
#'  and analyses of the joint occupancy decline and the parametric forms of this decline, for
#'   one particular community. In particular, \link[msco]{Jo.eng}:
#' \itemize{
#' \item{computes the joint occupancy (i.e. the number of sites or assemblages
#'  harbouring multiple species simultaneously)};
#'  \item{performs a null model test using the same index};
#'  \item{fits the three regression models (exponential, power law and exponential-power law)
#'   to joint occupancy decline (*sensu* Lagat *et al.,* 2021a) with order (number of species)};
#'   \item{estimates the parameter values of these models};
#'   \item{determines the best model among the three using AIC values};
#'   \item{quantifies the performance of the fitted models using the Pearson's \eqn{r^2}};
#'   \item{plots the joint occupancy decline regression and null models}, and
#'   \item{ascertains the archetypes of the patterns of species co-occurrences (from null
#'    model test) from which inferences on the type of drivers structuralising ecological
#'     communities can be made.}
#'}
#'
#' @param s.data A species-by-site presence/absence matrix with entries indicating
#' occurrence (1) and non-occurrence (0) of species in a site.
#' @param algo Randomisation algorithm used for the comparison with the null model. The
#'  possible options to choose from are: `sim1`, `sim2`, `sim3`, `sim4`, `sim5`, `sim6`,
#'   `sim7`, `sim8`, and `sim9`, all from Gotelli (2000). `sim2` is highly recommended
#'    (see Lagat *et al.,* 2021a).
#' @param metric The type of rescaling applied to the joint occupancy metric. Available options are:
#'  `Simpson_eqn` for Simpson equivalent, `Sorensen_eqn` for Sorensen equivalent, and `raw` for the
#'   raw form of index without rescaling.
#' @param nReps Number of simulations used in the null model test.
#' @param nmod_stats A Boolean indicating whether the summary
#'  statistics for the null model test should be output.
#' @param s.dplot A Boolean indicating whether the standard deviation
#'  of multi-species co-occurrence index should be included in the plots of joint occupancy
#'   decline or not.
#' @param All.plots A Boolean indicating whether joint occupancy decline  regression
#'  and null model plots should be output.
#' @param p.n.plot A Boolean indicating whether null model plot produced using the
#'  pairwise natural metric should be output.
#' @param trans A Boolean indicating if the observed and simulated values used in
#'  `p.n.plot` should be transformed by raising them to (1/100). This can be done to
#'   compare `p.n.plot` with `All.plots` at a point where the order, `i = 2`.
#' @param m.n.plot A Boolean indicating whether null model plot produced using joint
#'  occupancy metrics should be output. The default is `FALSE`.
#' @param dig The number of decimal places of the joint occupancy values (y axis) in the plots.
#'  The default is 3.
#' @param Jo.coeff A Boolean indicating if coefficient estimates of the joint occupancy
#'  decline regression models should be printed.
#' @param my.AIC A Boolean indicating whether Akaike Information Criterion of the joint occupancy
#'  decline regression models should be output or not.
#' @param my.rsq A Boolean indicating whether square of correlation coefficient between the
#'  observed and predicted values of joint occupancy should be output.
#' @param Exp_Reg A Boolean indicating if exponential regression parametric model should be
#'  printed.
#' @param P.law_Reg A Boolean indicating if power law regression parametric model should be printed.
#' @param Exp_p.l_Reg A Boolean indicating if exponential-power law regression parametric model
#'  should be printed.
#' @param Obs.data A Boolean indicating if observed/empirical data should be output.
#' @param Sim.data A Boolean indicating if simulated/random data produced using any of the
#'  simulation algorithms should be output.
#' @param Jo_val.sim A Boolean indicating if joint occupancy values of the
#'  simulated species-by-site presence/absence matrices should be output.
#' @param C.I_Jo_val.sim A Boolean indicating if 95% confidence interval of the joint occupancy values of the
#'  simulated data should be printed. This interval is the area under the null model.
#' @param Jo_val.obs A Boolean indicating if joint occupancy values of the
#'  observed species-by-site presence/absence matrices should be output.
#' @param Metric A Boolean indicating if metric used should be printed.
#' @param Algorithm A Boolean indicating if simulation algorithm used should be printed.
#' @param S.order A Boolean indicating if the number of species whose joint occupancy is computed
#'  should be printed.
#' @param Pt_Arch_Vals A Boolean indicating if character strings indicating the location of
#'  joint occupancy value of the observed data relative to the critical
#'   values of the 95% closed confidence interval for every order (number of species), should
#'    be printed.
#' @param Atype A Boolean indicating if a character string indicating the overall archetype of
#'  joint occupancy decline should be printed.This value must be \eqn{\in \{}"A1", "A2", "A3", "A4",
#'   "A5", "A6", "A7", "A8", "A9"\eqn{\}} or "NA". "NA" could be the combinations of two or more
#'    of the nine expected archetypes.
#' @param leg A Boolean indicating if the legend should be added to the `m.n.plot`. This parameter
#'  helps to control the appearance of plots in this function.
#' @param lab A Boolean indicating if the plot labels should be added to the `m.n.plot`. This parameter
#'  helps to control the appearance of plots in this function.
#'
#' @return `Jo.eng` function returns a list containing the following outputs:
#' \item{`all.plots`}{Joint occupancy decline regression and null model plots.}
#' \item{`jo.coeff`}{Coefficient estimates of the joint occupancy decline regression models.}
#' \item{`AIC`}{Akaike information criterion of the joint occupancy decline regression models.}
#' \item{`r2`}{Square of correlation coefficient between the observed and predicted values of
#'  joint occupancy.}
#' \item{`Exp_reg`}{Exponential regression parametric model.}
#' \item{`P.law_reg`}{Power law regression parametric model.}
#' \item{`Exp_p.l_reg`}{Exponential-power law regression parametric model.}
#' \item{`Obs.data`}{Observed/empirical data.}
#' \item{`Sim.data`}{Simulated/random data produced using any of the simulation algorithms.}
#' \item{`jo.val.sim`}{Joint occupancy value of the simulated species-by-site
#' presence/absence matrices.}
#' \item{`C.I_Jo_val.sim`}{95% confidence interval of the joint occupancy value of the simulated data.}
#' \item{`jo.val.obs`}{joint occupancy value of the observed species-by-site
#' presence/absence matrices.}
#' \item{`Metric`}{Metric used. It must be "`j.occ`".}
#' \item{`Algorithm`}{Simulation algorithm used.}
#' \item{`nReps`}{Number of simulations performed. This value together with the joint occupancy
#'  value of the observed data, constitutes the sampling distribution.}
#' \item{`s.order`}{Number of species whose joint occupancy is computed.}
#' \item{`Pt_Arch_vals`}{Character strings indicating the location of joint occupancy value
#'  of the observed data relative to the critical values of the 95% closed confidence interval
#' of the simulated data, for every order (number of species).}
#' \item{Archetype}{A character string indicating the overall archetype from `Pt_Arch_vals`.
#'  It must be \eqn{\in \{}"A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9"\eqn{\}} or
#'   "NA". "NA" could be the combinations of two or more
#'    of the nine expected archetypes (see \link[msco]{Arch_schem}).}
#' @references
#' \enumerate{
#'  \item{Lagat, V. K., Latombe, G. and Hui, C. (2021a). *A multi-species co-occurrence
#'  index to avoid type II errors in null model testing*. DOI: `<To be added>`.}
#'
#'  \item{Gotelli, N. J. (2000). Null model analysis of species co-occurrence patterns.
#'  *Ecology, 81(9)*, 2606-2621. <https://doi.org/10.1890/0012-9658(2000)081[2606:NMAOSC]2.0.CO;2>}
#' }
#' @examples
#' ex.data <- read.csv(system.file("extdata", "274.csv", package = "msco"))
#' j.en <- msco::Jo.eng(ex.data, algo="sim2", metric = "raw", nReps = 999,
#'            dig = 3, s.dplot = FALSE, All.plots = TRUE, Jo.coeff = TRUE,
#'            my.AIC = TRUE, my.rsq = TRUE, Exp_Reg = TRUE, P.law_Reg = TRUE,
#'            Exp_p.l_Reg = TRUE, Obs.data = FALSE, Sim.data = FALSE,
#'            Jo_val.sim = FALSE, C.I_Jo_val.sim = FALSE, Jo_val.obs = TRUE,
#'            Metric = TRUE, Algorithm = TRUE, S.order = TRUE,
#'            nmod_stats = TRUE, Pt_Arch_Vals = TRUE, Atype = TRUE,
#'            p.n.plot = FALSE, trans = FALSE, lab=FALSE, leg=FALSE, m.n.plot = FALSE)
#' j.en
#'
#' @export
#' @md

Jo.eng<-function(s.data, algo="sim2", metric = "raw", nReps = 999, dig = 3,
                     s.dplot = FALSE, All.plots = TRUE, Jo.coeff = TRUE, my.AIC = TRUE,
                     my.rsq = TRUE, Exp_Reg = TRUE, P.law_Reg = TRUE, Exp_p.l_Reg = TRUE,
                     Obs.data = FALSE, Sim.data = FALSE, Jo_val.sim = FALSE,
                     C.I_Jo_val.sim = FALSE, Jo_val.obs = TRUE, Metric = TRUE,
                     Algorithm = TRUE, S.order = TRUE, nmod_stats = TRUE, Pt_Arch_Vals = TRUE,
                     Atype = TRUE, p.n.plot = FALSE, trans = FALSE, lab=FALSE, leg=FALSE, m.n.plot = FALSE){

  if (!is.matrix(s.data)) {
    s.data <- as.matrix(s.data)
  }
  s.data <- s.data[rowSums(s.data) > 0, ] ## Remove rows with no species

  ######### Algorithms ##############

  sim1 <- function(s.data){
    matrix(sample(s.data), ncol = ncol(s.data))
  }
  sim2 <- function(s.data){
    t(apply(s.data, 1, sample))
  }
  sim3 <- function(s.data){
    apply(s.data, 2, sample)
  }
  vector_sample <- function(s.data, weights){
    x <- mat.or.vec(length(s.data), 1)
    x[sample(1:length(s.data), size = sum(s.data),
             prob = weights)] <- 1
    return(x)
  }
  sim4 <- function(s.data){
    t(apply(s.data, 1, vector_sample, weights = colSums(s.data)))
  }
  sim5 <- function(s.data){
    apply(s.data, 2, vector_sample, weights = rowSums(s.data))
  }
  sim6 <- function(s.data){
    matrixWeights <- outer(rep(1, nrow(s.data)), colSums(s.data))
    out <- matrix(vector_sample(s.data, weights = matrixWeights),
                  ncol = ncol(s.data))
  }
  sim7 <- function(s.data){
    matrixWeights <- outer(rowSums(s.data), rep(1, ncol(s.data)))
    matrix(vector_sample(s.data, weights = matrixWeights),
           ncol = ncol(s.data))
  }
  sim8 <- function(s.data){
    matrixWeights <- outer(rowSums(s.data), colSums(s.data))
    matrix(vector_sample(s.data, weights = matrixWeights),
           ncol = ncol(s.data))
  }
  sim9 <- function(s.data){
    sim9_single <- function(mydata){
      sam.rows <- sample.int(nrow(mydata), 2)
      mypair <- mydata[sam.rows, ]
      sum.1 = colSums(mypair) == 1
      if (sum(sum.1) > 1) {
        columns <- which(sum.1)
        mydata[sam.rows, columns] <- mypair[, sample(columns)]
      }
      return(mydata)
    }
    r.data <- s.data
    for (i in 1:1000) {
      r.data <- sim9_single(r.data)
    }
    return(r.data)
  }

  ############### Jo metric and its SD ######################

  jo.ex<-function(s.data, order, metric = metric){

    s.data <- as.matrix(s.data)
    richness <- colSums(s.data)
    p <- exp(lchoose(richness, order) - lchoose(nrow(s.data),order))
    jo.val <- sum(p)
    if(metric == "raw"){
      return(jo.val)
    }else if(metric=="Simpson_eqn"){
      simp.jo <- jo.val/min(rowSums(s.data))
      return(simp.jo)
    }else if(metric=="Sorensen_eqn"){
      sore.jo <- jo.val/mean(rowSums(s.data))
      return(sore.jo)
    }else if((metric %in% c("raw", "Simpson_eqn", "Sorensen_eqn"))!=TRUE){
      stop("Wrong option for the joint occupancy metric provided. It must either be 'raw',
           'Simpson_eqn', or 'Sorensen_eqn'.")
    }
  }
  jo.exps<-function(s.data, orders = 1:nrow(s.data), metric = metric){

    jo=0
    for (i in orders) {
      jo[i]=jo.ex(s.data, order = i, metric = metric)
    }
    return(jo)
  }

  jo.sd<-function(s.data, order){

    s.data <- as.matrix(s.data)
    richness <- colSums(s.data)
    similarity_mat <- t(s.data) %*% s.data
    p <- exp(lchoose(richness, order) - lchoose(nrow(s.data),order))
    covmat <- exp(lchoose(similarity_mat, order) - lchoose(nrow(s.data),order))

    for (j in 1:ncol(s.data)) {
      for (k in 1:ncol(s.data)) {
        covmat[j, k] <- covmat[j, k] - p[j] * p[k]
      }
    }
    jo.var <- choose(nrow(s.data), order)/
      (choose(nrow(s.data), order) - 1) * sum(covmat)

    jo.sd <- sqrt(jo.var)
    return(jo.sd)
  }

  jo.sds<-function(s.data, orders = 1:nrow(s.data)){

    SDs=0
    for (i in orders) {
      SDs[i]=jo.sd(s.data, order = i)
    }
    return(SDs)
  }

  if (suppressWarnings(is.na(as.numeric(s.data[2,1])))) {
    s.data <- s.data[, -1]
    class(s.data) <- "numeric"
  }

  #####################################################################
  ### Matrix of joint occupancy for different orders (as columns) #####
  ###                      from random matrices                   #####
  #####################################################################

  algoF <- get(algo)
  x <- lapply(1:nReps,
              function(x) algoF(s.data)) # Produce nReps random matrices
  x <- c(list(s.data), x)
  simulated.metric <- matrix(NA, nrow = (nReps+1), ncol = nrow(s.data))
  for (j in 1:nrow(s.data)) {
    for (i in 1:(nReps+1)) {
      sim.mat <- as.matrix(x[[i]])  # Select ith matrix from 'x' list
      simulated.metric[i,j] <- jo.ex(sim.mat, order=j, metric = metric)
    }
  }

  Sim <- simulated.metric  # joint occupancy (diff orders) of
  # nReps simulated data

  Obs <- jo.exps(s.data,1:nrow(s.data), metric = metric) # joint occupancy (diff orders)
  # of observed data

  SDs <- jo.sds(s.data, 1:nrow(s.data)) # s.d's of joint occupancy
  # (diff orders) of observed data

  ################# 95% Confidence Interval from Sim ###################

  C.I <- apply(Sim, 2, stats::quantile, probs=c(0,0.025,0.975,1), na.rm=TRUE)
  Min <- C.I[1,] # Min. value
  L.L <- C.I[2,] # Lower Critical Value
  U.L <- C.I[3,] # Upper Critical Value
  Max <- C.I[4,] # Max. value

  c.i <- apply(Sim, 2, stats::quantile, probs=c(0.025,0.975), na.rm=TRUE)

  ### Null model plot using pairwise natural metric  ####################

  if(p.n.plot == TRUE & trans==FALSE){
    grDevices::pdf(file = paste0(system.file("ms", package = "msco"), "/pairwise.nm.plot.pdf"), height = 5, width = 6)
    simm <- Sim[,2]
    obss <- Obs[2]
    ci <- c.i[,2]

    graphics::hist(simm, xlab = expression(paste("Simulated", " ",
                                                 italic(natural), " ",
                                                 "metric", " ", (italic(J)^{'{2}'}), sep='')),
                   main = "Null model test")
    graphics::abline(v=obss,col="blue", lty="dotted", lwd=1.8)
    graphics::abline(v=ci,col="red", lty="longdash", lwd=1.5)
    grDevices::dev.off()

  }else if(p.n.plot == TRUE & trans==TRUE){
    grDevices::pdf(file = paste0(system.file("ms", package = "msco"), "/pairwise.nm.plot.pdf"), height = 5, width = 6)
    simm <- round((Sim[,2])^(1/100), digits = dig)
    obss <- round((Obs[2])^(1/100), digits = dig)
    ci <-round((c.i[,2])^(1/100), digits = dig)
    graphics::hist(simm, xlab = expression(paste("Simulated", " ",
                                                 italic(natural), " ",
                                                 "metric", " ", (italic(J)^{'{2}'}), sep='')),
                   main = "Null model test")
    graphics::abline(v=obss,col="blue", lty="dotted", lwd=1.8)
    graphics::abline(v=ci,col="red", lty="longdash", lwd=1.5)
    grDevices::dev.off()
  }

  ###########################################################################
  #################### (e) joint occupancy Null model plot #####################
  ###########################################################################

  s.order <- 1:nrow(s.data)
  dtf <- data.frame(Obs, SDs, L.L, U.L)
  dtfl <- log10(dtf)
  ldtf <- data.frame(s.order, dtfl)
  ldtf <- do.call(data.frame, lapply(ldtf, function(x){
    replace(x, is.infinite(x) | is.na(x), NA)}))
  colnames(ldtf) <- c("s.order", "lObs", "lSDs", "lL.L", "lU.L")
  lObs <- ldtf$lObs
  ldtf <- ldtf[stats::complete.cases(ldtf),]
  lL.L <- ldtf$lL.L
  lU.L <- ldtf$lU.L

  nullmod.plot <- ggplot2::ggplot(ldtf, ggplot2::aes(x=s.order, y=lObs)) +
    ggplot2::theme_grey() +
    ggplot2::geom_point(ggplot2::aes(y = lObs, colour="Observed"),
                        shape = 1) +
    ggplot2::geom_line(ggplot2::aes(y = lObs, colour="Observed"),
                       size=0) +
    ggplot2::geom_line(ggplot2::aes(x = 2), colour="darkgrey",
                       size=0, linetype="solid") +
    ggplot2::geom_ribbon(data = ldtf, ggplot2::aes(ymin = lL.L,
                                                   ymax=lU.L,
                                                   fill="Null model"),
                         alpha=0.5) +
    ggplot2::scale_colour_manual(c("",""),values=c("black"))+
    ggplot2::scale_fill_manual("",values="grey12", drop=F)+
    ggplot2::ylab("")+
    ggplot2::xlab("") +
    ggplot2::scale_x_continuous(
      breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))+
    ggplot2::scale_y_continuous(labels = function(y)round((10^y)^(1/100),
                                                          digits = dig),
                                breaks = scales::pretty_breaks(n=4))+
    ggplot2::ggtitle("Null model test")+
    ggplot2::theme(legend.position = "none") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  ### Archetypes

  Archetype <- dplyr::case_when(
    ldtf$lObs > ldtf$lU.L ~ "A1",
    ldtf$lObs < ldtf$lL.L ~ "A9",
    ldtf$lObs <= ldtf$lU.L & ldtf$lObs >= ldtf$lL.L ~ "A5"
  )

  myArch <- c()
  if(all(Archetype[2:length(Archetype)]=="A1")==TRUE){
    myArch <- "A1"
  }else if(all(all(Archetype[which(Archetype=="A1")]=="A1")==TRUE &
               length(Archetype[which(Archetype=="A1")]=="A1")>=1 &
               all(Archetype[(length(which(
                 Archetype=="A1"))+2):length(Archetype)]=="A5")==TRUE)==TRUE){
    myArch = "A2"
  }else if(all(all(Archetype[which(Archetype=="A1")]=="A1")==TRUE &
               length(Archetype[which(Archetype=="A1")]=="A1")>=1 &
               all(Archetype[(length(which(
                 Archetype=="A1"))+2):length(Archetype)][which(Archetype[(length(which(
                   Archetype=="A1"))+2):length(Archetype)]=="A5")]=="A5")==TRUE &
               length(Archetype[(length(which(
                 Archetype=="A1"))+2):length(Archetype)][which(Archetype[(length(which(
                   Archetype=="A1"))+2):length(Archetype)]=="A5")]=="A5")>=1 &
               all(Archetype[(length(which(
                 Archetype=="A1"))+2):length(Archetype)][(length(which(
                   Archetype=="A5"))+1):length(Archetype[(length(which(
                     Archetype=="A1"))+2):length(Archetype)])]=="A9")==TRUE)==TRUE){
    myArch = "A3"
  }else if(all(all(Archetype[which(Archetype=="A5")]=="A5")==TRUE &
               length(Archetype[which(Archetype=="A5")]=="A5")>=2 &
               all(Archetype[(length(which(
                 Archetype=="A5"))+1):length(Archetype)]=="A1")==TRUE)==TRUE){
    myArch = "A4"
  }else if(all(Archetype[1:length(Archetype)]=="A5")==TRUE){
    myArch = "A5"
  }else if(all(all(Archetype[which(Archetype=="A5")]=="A5")==TRUE &
               length(Archetype[which(Archetype=="A5")]=="A5")>=2 &
               all(Archetype[(length(which(
                 Archetype=="A5"))+1):length(Archetype)]=="A9")==TRUE)==TRUE){
    myArch = "A6"
  }else if(all(all(Archetype[which(Archetype=="A9")]=="A9")==TRUE &
               length(Archetype[which(Archetype=="A9")]=="A9")>=1 &
               all(Archetype[(length(which(
                 Archetype=="A9"))+2):length(Archetype)][which(Archetype[(length(which(
                   Archetype=="A9"))+2):length(Archetype)]=="A5")]=="A5")==TRUE &
               length(Archetype[(length(which(
                 Archetype=="A9"))+2):length(Archetype)][which(Archetype[(length(which(
                   Archetype=="A9"))+2):length(Archetype)]=="A5")]=="A5")>=1 &
               all(Archetype[(length(which(
                 Archetype=="A9"))+2):length(Archetype)][(length(which(
                   Archetype=="A5"))+1):length(Archetype[(length(which(
                     Archetype=="A9"))+2):length(Archetype)])]=="A1")==TRUE)==TRUE){
    myArch = "A7"
  }else if(all(all(Archetype[which(Archetype=="A9")]=="A9")==TRUE &
               length(Archetype[which(Archetype=="A9")]=="A9")>=1 &
               all(Archetype[(length(which(
                 Archetype=="A9"))+2):length(Archetype)]=="A5")==TRUE)==TRUE){
    myArch = "A8"
  }else if(all(Archetype[2:length(Archetype)]=="A9")==TRUE){
    myArch = "A9"
  }else{
    myArch = "NA"
  }


  ###########################################################################
  ##################### joint occupancy decline Regression plots ###############
  ###########################################################################

  ##################### (a) joint occupancy decline ############################

  lObs <- log10(Obs)
  jod <- data.frame(s.order,lObs)
  jod <- do.call(data.frame, lapply(jod, function(x){
    replace(x, is.infinite(x) | is.na(x), NA)}))
  jod <- jod[stats::complete.cases(jod),]
  # lObs <- jod$lObs
  # s.order <- jod$s.order
  p1 <- ggplot2::ggplot(jod, ggplot2::aes(x=s.order, y=lObs))+
    ggplot2::theme_grey() +
    ggplot2::geom_point(ggplot2::aes(y = lObs), shape = 1) +
    ggplot2::xlab("")+
    ggplot2::ylab("") +
    ggplot2::scale_x_continuous(
      breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))+
    ggplot2::scale_y_continuous(labels = function(y)round((10^y)^(1/100),
                                                          digits = dig),
                                breaks = scales::pretty_breaks(n=4))+
    ggplot2::ggtitle("J. occ. decl.")+
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::theme(legend.position = "none")


  if(s.dplot == TRUE){

    upper.sd <- Obs + SDs
    lower.sd <- Obs - SDs
    dtf.sd <- data.frame(Obs, upper.sd, lower.sd)
    dtfl.sd <- dtf.sd
    ldtf.sd <- data.frame(s.order, dtfl.sd)
    ldtf.sd <- do.call(data.frame, lapply(ldtf.sd, function(x){
      replace(x, is.infinite(x) | is.na(x), NA)}))
    colnames(ldtf.sd) <- c("s.order.sd", "lObs", "lupper.sd", "llower.sd")

    ldtf.sd <- ldtf.sd[stats::complete.cases(ldtf.sd),]

    row_sub = apply(ldtf.sd, 1, function(row) all(row !=0 ))
    ldtf.sd <- ldtf.sd[row_sub,]

    lObs.sd <- ldtf.sd$lObs
    lupper.sd <- ldtf.sd$lupper.sd
    llower.sd <- ldtf.sd$llower.sd
    s.order.sd <- ldtf.sd$s.order.sd

    p1 <- ggplot2::ggplot(ldtf.sd, ggplot2::aes(x=s.order.sd, y=lObs.sd))+
      ggplot2::theme_grey() +
      ggplot2::geom_point(ggplot2::aes(y = lObs.sd), shape = 1) +
      ggplot2::geom_line(ggplot2::aes(y = lupper.sd), size=0.1) +
      ggplot2::geom_line(ggplot2::aes(y = llower.sd), size=0.1) +

      ggplot2::scale_colour_manual("",values=c("black"))+
      ggplot2::xlab("")+
      ggplot2::ylab("") +
      ggplot2::scale_x_continuous(
        breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))+
      ggplot2::scale_y_continuous(labels = function(y)round((10^y)^(1/100),
                                                            digits = dig),
                                  breaks = scales::pretty_breaks(n=4))+
      ggplot2::ggtitle("J. occ. decl.")+
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::theme(legend.position = "none")
  }


  ##################### (b) Exponential regression #############################

  fm <- stats::lm(jod$lObs ~ jod$s.order, na.action = stats::na.omit)
  jo.exp <- suppressWarnings(minpack.lm::nlsLM(Obs ~ p * exp(q * s.order),
                                               start =list(p=exp(stats::coef(fm)[1]), q=stats::coef(fm)[2])))
  pred.exp <- log10(stats::predict(jo.exp, data.frame(s.order)))
  pred.exp[is.infinite(pred.exp)]<-NA
  lObs <- log10(Obs)
  lObs[is.infinite(lObs)]<-NA
  jode <- data.frame(s.order, lObs, pred.exp)
  jode <- jode[stats::complete.cases(jode),]

  p2 <- ggplot2::ggplot(jode, ggplot2::aes(x=s.order, y=lObs)) +
    ggplot2::theme_grey() +
    ggplot2::geom_point(ggplot2::aes(y = lObs, colour="Observed"), shape = 1) +
    ggplot2::geom_line(ggplot2::aes(y = pred.exp, colour="Predicted"),
                       linetype = "dashed", size=1) +
    ggplot2::scale_colour_manual("",values=c("black", "black"))+
    ggplot2::ylab("") +
    ggplot2::xlab("") +
    ggplot2::scale_x_continuous(
      breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))+
    ggplot2::scale_y_continuous(labels = function(y)round((10^y)^(1/100),
                                                          digits = dig),
                                breaks = scales::pretty_breaks(n=4))+
    ggplot2::ggtitle("Exp. reg.")+
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::theme(legend.position = "none")


  ##################### (c) Power-law regression ##############################

  ls.order = log10(jod$s.order)
  fmp <- stats::lm(jod$lObs ~ ls.order, na.action = stats::na.omit)
  jo.pl <- suppressWarnings(minpack.lm::nlsLM(Obs ~ m * (s.order)^n,
                                              start =list(m=1, n=stats::coef(fmp)[2])))
  pred.pl <- log10(stats::predict(jo.pl, data.frame(s.order)))
  pred.pl[is.infinite(pred.pl)]<-NA
  lObs <- log10(Obs)
  lObs[is.infinite(lObs)]<-NA
  jodpl <- data.frame(s.order,lObs, pred.pl)
  jodpl <- jodpl[stats::complete.cases(jodpl),]

  p3 <- ggplot2::ggplot(jodpl, ggplot2::aes(x=s.order, y=lObs)) +
    ggplot2::theme_grey() +
    ggplot2::geom_point(ggplot2::aes(y = lObs, colour="Observed"),
                        shape = 1) +
    ggplot2::geom_line(ggplot2::aes(y = pred.pl, colour="Predicted"),
                       linetype = "dashed", size = 1) +
    ggplot2::scale_colour_manual("",values=c("black", "black"))+
    ggplot2::xlab("") +
    ggplot2::ylab("") +
    ggplot2::scale_x_continuous(
      breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))+
    ggplot2::scale_y_continuous(labels = function(y)round((10^y)^(1/100),
                                                          digits = dig),
                                breaks = scales::pretty_breaks(n=4))+
    ggplot2::ggtitle("P.law reg.")+
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::theme(legend.position = "none")



  ##################### (d) Exponential-Power law regression #####################


  jo.exp.pl <- suppressWarnings(minpack.lm::nlsLM(Obs~a*exp(b*s.order)*s.order^c,
                                                  start=list(a=1, b=0, c=0)))
  pred.exp.pl <- log10(stats::predict(jo.exp.pl, data.frame(s.order)))
  pred.exp.pl[is.infinite(pred.exp.pl)]<-NA
  lObs <- log10(Obs)
  lObs[is.infinite(lObs)]<-NA
  jodepl <- data.frame(s.order,lObs, pred.exp.pl)
  jodepl <- jodepl[stats::complete.cases(jodepl),]

  p4 <- ggplot2::ggplot(jodepl, ggplot2::aes(x=s.order, y=lObs)) +
    ggplot2::theme_grey() +
    ggplot2::geom_point(ggplot2::aes(y = lObs, colour="Observed"), shape = 1) +
    ggplot2::geom_line(ggplot2::aes(y = pred.exp.pl, colour="Predicted"),
                       linetype = "dashed", size=1) +
    ggplot2::scale_colour_manual("",values=c("black", "black"))+
    ggplot2::xlab("") +
    ggplot2::ylab("") +
    ggplot2::scale_x_continuous(
      breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))+
    ggplot2::scale_y_continuous(labels = function(y)round((10^y)^(1/100),
                                                          digits = dig),
                                breaks = scales::pretty_breaks(n=4))+
    ggplot2::ggtitle("E-P.l. reg.")+
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))+
    ggplot2::theme(legend.position = "none")

  ###############################################################################
  ##########                                                              #######
  ########## Parameter values & goodness of fit for all model regressions #######
  ##########                                                              #######
  ###############################################################################

  trunc <- function(x, ..., prec = 0) base::trunc(x * 10^prec, ...) / 10^prec
  rsq <- function (x, y) stats::cor(x, y) ^ 2

  pred.exp <- stats::predict(jo.exp, data.frame(s.order))
  pred.pl <- stats::predict(jo.pl, data.frame(s.order))
  pred.exp.pl <- stats::predict(jo.exp.pl, data.frame(s.order))

  exp.r2 <- trunc(rsq(Obs, pred.exp),prec=4)
  pl.r2 <-  trunc(rsq(Obs, pred.pl),prec=4)
  exp.pl.r2 <- trunc(rsq(Obs, pred.exp.pl),prec=4)

  rsq <- c(exp.r2, pl.r2, exp.pl.r2)

  ########## AIC for exp, p.law & exp-p.law models ##################

  aic1 <- stats::AIC(jo.exp, jo.pl, jo.exp.pl)
  Delta_AIC3 <- (aic1[,2])-min(aic1[,2])

  aic3 <- stats::AIC(jo.exp, jo.pl)
  Delta_AIC2 <- (aic3[,2])-min(aic3[,2])

  Delta_AIC2 <- c(Delta_AIC2, NA)
  aic <- cbind(aic1, Delta_AIC3, Delta_AIC2)




  jo.exp.coeff <-data.frame(t(as.matrix(c(stats::coef(jo.exp)[[1]],
                                          stats::coef(jo.exp)[[2]], NA))))

  jo.pl.coeff <- data.frame(t(as.matrix(c(stats::coef(jo.pl)[[1]],
                                          stats::coef(jo.pl)[[2]], NA))))

  jo.exp.pl.coeff <- data.frame(t(as.matrix(c(stats::coef(jo.exp.pl)[[1]],
                                              stats::coef(jo.exp.pl)[[2]],
                                              stats::coef(jo.exp.pl)[[3]]))))

  jo.coeff <- as.data.frame(rbind(jo.exp.coeff,jo.pl.coeff, jo.exp.pl.coeff))

  rownames(jo.coeff) <- c("Exponential", "Power law", "Exp-p. law")
  colnames(jo.coeff) <- c("a", "b", "c")

  #####################################################################
  ######################### All plots #################################
  #####################################################################

  ex.aic <- ceiling(aic[1,2])
  pl.aic <- ceiling(aic[2,2])
  ex.pl.aic <- ceiling(aic[3,2])

  all.plots <- cowplot::ggdraw() +
    cowplot::draw_plot(p1, x = 0.07, y = .55, width = .22, height = .45) +
    cowplot::draw_plot(p2, x = .29, y = 0.55, width = .22, height = .45) +
    cowplot::draw_plot(p3, x = 0.07, y = 0.1, width = .22, height = .45) +
    cowplot::draw_plot(p4, x = 0.29, y = 0.1, width = .22, height = .45) +
    cowplot::draw_plot(nullmod.plot, x = .51, y = 0.1, width = 0.36,
                       height = 0.9) +

    cowplot::draw_plot_label(label = c("(a)", "(b)", "(c)", "(d)", "(e)", myArch),
                             size = 12, x = c(0.11, 0.34, 0.12, 0.34, 0.56,
                                              0.8), y = c(1, 1, 0.55, 0.55,
                                                          1, 0.92)) +

    cowplot::draw_plot_label(label = c(paste("AIC =",ex.aic),
                                       paste("AIC =",pl.aic),
                                       paste("AIC =",ex.pl.aic)),
                             size = 8, x = c(0.438, 0.222, 0.438),
                             y = c(0.95, 0.495, 0.495),hjust = 0)+

    cowplot::draw_plot_label(label = c(paste("rsq >",trunc(100 * exp.r2) / 100),
                                       paste("rsq >",trunc(100 * pl.r2) / 100),
                                       paste("rsq >",trunc(100 * exp.pl.r2) / 100)),
                             size = 8, x = c(0.438, 0.222, 0.438),
                             y = c(0.925, 0.47, 0.47),hjust = 0)+

    cowplot::draw_image(system.file("logos", "key.png", package = "msco"),
                        x=0.42, y=0.07,scale = 0.2) +

    cowplot::draw_image(system.file("logos", "ylab.jpg", package = "msco"),
                        x=0.02, y=0.09, scale = 0.65,
                        width = 0.05) +

    cowplot::draw_image(system.file("logos", "xlab.jpg", package = "msco"),
                        x=0.0, y=0.02, scale = 0.57,
                        height = 0.07)


  ##### null model plot only


  if(leg==TRUE){
    nplot1 <- ggplot2::ggplot(ldtf, ggplot2::aes(x=s.order,
                                                 y=lObs)) +
      ggplot2::theme_grey() +
      ggplot2::geom_point(ggplot2::aes(y = lObs, colour="Observed"),
                          shape = 1) +
      ggplot2::geom_line(ggplot2::aes(y = lObs, colour=c("Observed")),
                         size=0) +
      ggplot2::geom_line(ggplot2::aes(x = 2), colour="darkgrey",
                         size=0, linetype="dashed") +
      ggplot2::geom_ribbon(data = ldtf, ggplot2::aes(ymin = lL.L,
                                                     ymax=lU.L,
                                                     fill="Null model"),
                           alpha=0.5) +
      ggplot2::scale_colour_manual(c("",""),values=c("black"))+
      ggplot2::scale_fill_manual("",values="grey12", drop=F)+
      ggplot2::ylab("")+
      ggplot2::xlab("") +
      ggplot2::scale_x_continuous(
        breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))+
      ggplot2::scale_y_continuous(labels = function(y)round((10^y)^(1/100),
                                                            digits = dig),
                                  breaks = scales::pretty_breaks(n=4))+
      ggplot2::ggtitle("")+
      ggplot2::theme(legend.position = "right")
  }else{
    nplot1 <- ggplot2::ggplot(ldtf, ggplot2::aes(x=s.order,
                                                 y=lObs)) +
      ggplot2::theme_grey() +
      ggplot2::geom_point(ggplot2::aes(y = lObs, colour="Observed"),
                          shape = 1) +
      ggplot2::geom_line(ggplot2::aes(y = lObs, colour=c("Observed")),
                         size=0) +
      ggplot2::geom_line(ggplot2::aes(x = 2), colour="darkgrey",
                         size=0, linetype="dashed") +
      ggplot2::geom_ribbon(data = ldtf, ggplot2::aes(ymin = lL.L,
                                                     ymax=lU.L,
                                                     fill="Null model"),
                           alpha=0.5) +
      ggplot2::scale_colour_manual(c("",""),values=c("black"))+
      ggplot2::scale_fill_manual("",values="grey12", drop=F)+
      ggplot2::ylab("")+
      ggplot2::xlab("") +
      ggplot2::scale_x_continuous(
        breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))+
      ggplot2::scale_y_continuous(labels = function(y)round((10^y)^(1/100),
                                                            digits = dig),
                                  breaks = scales::pretty_breaks(n=4))+
      ggplot2::ggtitle("")+
      ggplot2::theme(legend.position = "none")
  }


  nplot <- cowplot::ggdraw() +
    if(lab==TRUE & leg==TRUE){
      cowplot::draw_plot(nplot1, x = 0.1, y = 0.1, width = 0.9, height = 0.9)+
      cowplot::draw_image(system.file("logos", "ylab.jpg", package = "msco"),
                            x=0.05, y=0.1, scale = 0.65, width = 0.05) +
      cowplot::draw_image(system.file("logos", "xlab.jpg", package = "msco"),
                            x=0.05, y=0.05, scale = 0.6, height = 0.05)+
      cowplot::draw_plot_label(label = myArch, size = 12, x=0.8, y = 0.93)
    }else if(lab==TRUE & leg==FALSE){
      cowplot::draw_plot(nplot1, x = 0.1, y = 0.1, width = 0.9, height = 0.9)+
      cowplot::draw_image(system.file("logos", "ylab.jpg", package = "msco"),
                            x=0.05, y=0.1, scale = 0.65, width = 0.05) +
      cowplot::draw_image(system.file("logos", "xlab.jpg", package = "msco"),
                            x=0.05, y=0.05, scale = 0.6, height = 0.05)
    }else if(lab==FALSE & leg==TRUE){
      cowplot::draw_plot(nplot1, x = 0.1, y = 0.1, width = 0.9, height = 0.9)
    }else if(lab==FALSE & leg==FALSE){
      cowplot::draw_plot(nplot1, x = 0.1, y = 0.1, width = 0.9, height = 0.9)
    }

  #######################################################################


  jo.engine <- list()

  if(All.plots == TRUE){
    jo.engine$all.plots <- all.plots
  }
  if(m.n.plot == TRUE){
    jo.engine$m.n.plot <- nplot
  }

  ## Pairwise null model plot
  if(p.n.plot==TRUE){
    jo.engine$p.n.plot <- print(noquote("Check msco's 'ms' folder in your R version's directory for a 'pairwise.nm.plot.pdf' file."))
  }

  jo.engine$"#################  REGRESSION ANALYSES  ###############" <- noquote('')

  if(Jo.coeff == TRUE){
    jo.engine$jo.coeff <- jo.coeff
  }
  if(my.AIC == TRUE){
    jo.engine$AIC <- aic
  }
  if(my.rsq == TRUE){
    jo.engine$r2 <- rsq
  }
  if(Exp_Reg == TRUE){
    jo.engine$Exp_reg <- jo.exp
  }
  if(P.law_Reg == TRUE){
    jo.engine$P.law_reg <- jo.pl
  }
  if(Exp_p.l_Reg == TRUE){
    jo.engine$Exp_p.l_reg <- jo.exp.pl
  }

  jo.engine$"#############  NULL MODEL ANALYSIS  #############" <- noquote("")

  if(Obs.data == TRUE){
    jo.engine$Obs.Data <- s.data
  }
  if(Sim.data == TRUE){
    jo.engine$Sim.Data <- sim.mat
  }
  if(Jo_val.sim == TRUE){
    jo.engine$jo.val.sim <- Sim
  }
  if(C.I_Jo_val.sim == TRUE){
    jo.engine$C.I_Jo_val.sim <- c.i
  }
  if(Jo_val.obs == TRUE){
    jo.engine$jo.val.obs <- Obs
  }
  if(Metric == FALSE){
    jo.engine$Index <- metric
  }
  if(Algorithm == FALSE){
    jo.engine$algorithm <- algo
  }
  if(S.order == TRUE){
    jo.engine$s.order <- s.order
  }

  ################ Nullmod Stats ###############
  Order <- s.order
  nmodd <- as.data.frame(rbind(Order, Obs, L.L, U.L))
  nmod <- `row.names<-`(nmodd, c('Order', 'Obs', 'L.C.V', 'U.C.V'))
  if(nmod_stats==TRUE){
    jo.engine$nmod_stats <- `colnames<-`(round(nmod,5), NULL)
  }

  if(Pt_Arch_Vals == TRUE){
    jo.engine$Pt_Arch_vals <- Archetype
  }
  if(Atype == TRUE){
    jo.engine$Archetype <- myArch
  }

  class(jo.engine) <- "jo.Obj"

  return(jo.engine)

}

#' A generalised B-spline modelling for a set of neutral and trait-based variables
#'
#' This function implements the generalised B-spline model (*sensu* Lagat *et al.,* 2021b)
#'  for dissecting the effects of random encounter versus functional trait mismatching on
#'   multi-species co-occurrence and interference. Generalized linear model
#'   (*sensu* Hastie and Tibshirani, 1986) with binomial variance distribution and log link
#'    functions employed, with predictors transformed using a linear combination of B-splines
#'     (*sensu* Curry and Schoenberg, 1988).
#'
#' @param s.data A species-by-site presence/absence `data.frame` with entries indicating
#' occurrence (1) and non-occurrence (0) of species in a site.
#' @param t.data A `data.frame` with traits as columns and species as rows. The species must be the same as in `s.data`.
#' @param p.d.mat A symmetric `matrix` with dimension names as species and entries indicating the
#'  phylogenetic distance between any two of them (species).
#' @param d.f Degrees of freedom for B-splines.
#' @param order.jo Specific number of species for which the joint occupancy is computed. To implement
#'  generalised B-spline modelling for multiple orders, see \link[msco]{gbsm_m.orders} function.
#' @param degree Degree of the B-splines.
#' @param metric The type of rescaling applied to the joint occupancy metric. Available options are:
#'  `Simpson_eqn` for Simpson equivalent, `Sorensen_eqn` for Sorensen equivalent, `raw_prop` for the
#'   raw form of the metric rescaled by dividing by the total number of sites, N, and `raw` for the
#'    raw form of the metric without rescaling.
#' @param gbsm.model The model used if the `raw` form of the metric is choosen. Availbale options are `"quasipoisson"`
#'   for quasipoisson GLM or `"nb"` for negative binomial GLM. Other metric types strictly uses binomial GLM.
#' @param n Number of samples for which the joint occupancy is computed. These samples are non-overlapping.
#'  I.e., sampling is done without replacement. If the total number of combinations of `i` species chosen
#'   from the total species pool `m`, i.e. `choose(m,i)`, is less than this value (`n`), `choose(m,i)` is
#'    used as the (maximum) number of samples one can set. Otherwise sampling without replacement is
#'     performed to select just the `n` samples.
#' @param b.plots Boolean value indicating if B-spline basis functions should be plotted.
#' @param bsplines This parameter indicates if a single or all B-spline curves should be plotted.
#'  If `b.plots=TRUE` and `bsplines="single"`, the B-splines for the first predictor (labelled, ` "Trait variable" `)
#'   in `t.data` will be plotted. Any other value for `bsplines` (other than `"single"`) results in the
#'    B-splines for all predictors being plotted.
#' @param response.curves A boolean value indicating if all response curves should be plotted.
#' @param leg Boolean value indicating if the legend of the gbsm outputs should be included in the plots. This
#'  parameter is added to help control the appearance of plots in \link[msco]{gbsm_m.orders} function.
#' @param ylabel Boolean value indicating if the y label should be included in the response curves. This
#'  parameter is added to help control the appearance of plots in \link[msco]{gbsm_m.orders} function.
#' @param scat.plot Boolean value indicating if scatter plots between joint occupancy and its predicted
#'  values should be plotted.
#' @param start.range Range of starting values for glm regression.
#'
#' @return `gbsm` function returns a list containing the following outputs:
#' \item{`order.jo`}{Order of joint occupancy}
#' \item{`Predictors`}{Predictor variables used in GLM regression with binomial variance
#'  distribution function and log link function.}
#' \item{`Responses`}{Response variables from GLM regression with binomial variance distribution
#'  function and log link function.}
#' \item{`coeff`}{Coefficients of the generalized linear model used.}
#' \item{`glm_obj`}{Generalized linear model used.}
#' \item{`j.occs`}{Rescaled observed joint occupancies. See `metric`above.}
#' \item{`bs_pred`}{B-spline-transformed `Predictors`.}
#' \item{`var.expld`}{Amount of variation in joint occupancy explained by the `Predictors`. I.e.,
#'  it is the Pearson's **\eqn{r^2}** between the observed and predicted values of joint occupancy.}
#' \item{`summary`}{summary of the regression results}
#'
#' @references
#' \enumerate{
#'
#'  \item{Curry, H. B., and Schoenberg, I. J. (1988). On P&oacute;lya frequency functions IV: the
#'   fundamental spline functions and their limits. In *IJ Schoenberg Selected Papers*
#'    (pp. 347-383). Birkh&auml;user, Boston, MA. <https://doi.org/10.1007/978-1-4899-0433-1_17>}
#'
#'  \item{Hastie, T., and Tibshirani, R. (1986). Generalized additive models. *Stat. Sci. 1*(3),
#'   297-310. <https://doi.org/10.1214/ss/1177013604>}
#'
#'  \item{Lagat, V. K., Latombe, G. and Hui, C. (2021a). *A multi-species co-occurrence
#'  index to avoid type II errors in null model testing*. DOI: `<To be added>`.}
#'
#'  \item{Lagat, V. K., Latombe, G. and Hui, C. (2021b). *Dissecting the effects of random
#'  encounter versus functional trait mismatching on multi-species co-occurrence and
#'   interference with generalised B-spline modelling*. DOI: `<To be added>`.}
#'  }
#'
#' @examples
#' \dontrun{
#'  my.path <- system.file("extdata/gsmdat", package = "msco")
#'  setwd(my.path)
#'  s.data <- get(load("s.data.csv")) ## Species-by-site matrix
#'  t.data <- get(load("t.data.csv")) ## Species-by-trait matrix
#'  p.d.mat <- get(load("p.d.mat.csv")) ## Species-by-species phylogenetic distance matrix
#'
#'  my.gbsm <- msco::gbsm(s.data, t.data, p.d.mat, metric = "Simpson_eqn", gbsm.model,
#'   d.f=4, order.jo=3, degree=3, n=1000, b.plots=TRUE, scat.plot=TRUE,
#'    bsplines="single", response.curves=TRUE, leg=1, start.range=c(-0.1,0))
#'
#'  my.gbsm$bs_pred
#'  my.gbsm$Predictors
#'  my.gbsm$Responses
#'  my.gbsm$order.jo
#'  my.gbsm$var.expld
#'  }
#'
#' @export
#' @md

gbsm <- function(s.data, t.data, p.d.mat, metric= "Simpson_eqn", d.f=4, order.jo=3, degree=3, n=1000, b.plots=TRUE,
                 bsplines="single", gbsm.model, scat.plot=TRUE, response.curves=TRUE, ylabel=TRUE, leg=1, start.range=c(-0.1,0)){

  if(class(t.data)!="data.frame"){
    t.data <- as.data.frame(t.data)
  }
  if (order.jo > dim(s.data)[1]) {
    stop("Wrong value for \"order\": it must be equal or lower than the number of species.")
  }
  if (!isSymmetric(as.matrix(p.d.mat))) {
    stop("Distance matrix is not symmetrical")
  }
  if(nrow(s.data)!=nrow(t.data)){
    stop("s.data and t.data must have the same number of rows")
  }
  if(nrow(s.data)!=nrow(p.d.mat)){
    stop("s.data and p.d.mat must have the same number of rows")
  }

  for (i in 1:ncol(t.data)) {
    t.data[, i] <- (t.data[, i] - min(t.data[, i]))/(max(t.data[, i])-min(t.data[, i]))
  }

  #### Compute the B-splines of the original trait variables (t.data) to get bs.traits
  bs.traits <- matrix(NA, nrow(t.data), (ncol(t.data) * d.f))
  for (j in 1:ncol(t.data)) {
    for (i in 1:d.f) {
      bs.traits[,(j - 1) * d.f + i] <- splines2::bSpline(t.data[,j], d.f=d.f, degree=degree, intercept = TRUE)[,i]
    }
  }
  bs.traits <- data.frame(bs.traits)

  ## Assign names to bs.traits
  for (i in 1:(ncol(bs.traits)/d.f)) {
    for (j in 1:d.f) {
      names(bs.traits)[(i - 1) * d.f + j] <- paste(names(t.data)[i], j, sep = "")
    }
  }
  bs.traits <- `rownames<-`(bs.traits, rownames(t.data))

  ## Compute the differences on the values transformed using b-splines (bs.traits) using SD, to get bs.traits.diff. Simultaneously compute p.dist and encounter rate
  order <- 1:nrow(s.data) ## Possible joint occupancy orders
  sn <- dim(s.data)[1] ## Total number of species
  N <- dim(s.data)[2] ## Total number of sites
  ncom <- choose(sn,order.jo) ## Total number of the combinations of "order.jo" species chosen from sn

  if(ncom > n){ ## Use MCMC to sample n species combinations if the total combinations > the chosen sample size, n.
    jo <- rep(NA, n)
    bs.traits.diff <- matrix(NA, nrow = n, ncol = ncol(bs.traits))
    erate <- rep(NA, n)
    p.dist <- rep(NA, n)

    for (j in 1:n) {
      sam <- sample(1:sn, order.jo, replace = FALSE)
      bs.traits.diff[j,] <- apply(bs.traits[sam,], 2, stats::sd)

      #### jo values for chosen combination of species
      jo[j] <- (msco::j.occ(s.data[sam,], order = order.jo)$jo.val)
      if(metric=="raw"){
        jo[j] <- jo[j]
      }else if(metric=="Simpson_eqn"){
        jo[j] <- jo[j]/min(rowSums(s.data[sam,]))
      }else if(metric=="Sorensen_eqn"){
        jo[j] <- jo[j]/mean(rowSums(s.data[sam,]))
      }
      ## p.dist
      p.dist[j] <- mean(p.d.mat[t(utils::combn(sort(sam), 2))])

      ## Encounter rate
      nsam <- length(sam)
      er <- rep(NA, nsam)
      for (i in 1:nsam) {
        er[i] <- rowSums(s.data[sam[i],])
      }
      erate[j] <- (prod(er))/(N^nsam)
    }
  }else{ ## Otherwise use the combinations as they are without sampling n of them.
    jo <- rep(NA, ncom)
    bs.traits.diff <- matrix(NA, nrow = ncom, ncol = ncol(bs.traits))
    erate <- rep(NA, ncom)
    p.dist <- rep(NA, ncom)
    com <- utils::combn(sn, order.jo)
    for (j in 1:ncom) {
      bs.traits.diff[j,] <- apply(bs.traits[com[,j],], 2, stats::sd)

      #### jo values for chosen combination of species
      if(metric=="raw"){
        jo[j] <- (msco::j.occ(s.data[com[,j],], order = order.jo)$jo.val)
      }else if(metric=="Simpson_eqn"){
        jo[j] <- (msco::j.occ(s.data[com[,j],], order = order.jo)$jo.val)/min(rowSums(s.data[com[,j],]))
      }else if(metric=="Sorensen_eqn"){
        jo[j] <- (msco::j.occ(s.data[com[,j],], order = order.jo)$jo.val)/mean(rowSums(s.data[com[,j],]))
      }

      ## p.dist
      p.dist[j] <- mean(p.d.mat[t(utils::combn(sort(com[,j]), 2))])

      ## Encounter rate
      nsam <- length(com[,j])
      er <- rep(NA, nsam)
      for (i in 1:nsam) {
        er[i] <- rowSums(s.data[com[,j][i],])
      }
      erate[j] <- (prod(er))/(N^nsam)
    }
  }

  if(metric=="raw_prop"){
    jo <- jo/N
  }
  if((metric %in% c("raw", "raw_prop", "Simpson_eqn", "Sorensen_eqn"))!=TRUE){
    stop("Wrong option for the joint occupancy metric provided. It must either be 'raw', 'raw_prop', 'Simpson_eqn', or 'Sorensen_eqn'.")
  }
  ## Rescale p.dist to be in [0,1]
  p.dist<- (p.dist - min(p.dist))/(max(p.dist)-min(p.dist))

  ## Assign bs.traits.diff (trait differences of B-splines of t.data) column names
  bs.traits.diff <- `colnames<-`(bs.traits.diff, colnames(bs.traits))
  bs.traits.diff <- data.frame(bs.traits.diff)

  ## B-splines of p.dist & erate
  bs.p.dist <- matrix(NA, nrow=n, ncol=d.f)
  bs.erate <- matrix(NA, nrow=length(erate), ncol=d.f)
  for (i in 1:d.f) {
    bs.p.dist[,i] <- splines2::bSpline(p.dist, d.f=d.f, degree=degree, intercept = TRUE)[,i]
    bs.erate[,i] <- splines2::bSpline(erate, d.f=d.f, degree=degree, intercept = TRUE)[,i]
  }
  bs.p.dist <- data.frame(bs.p.dist)
  bs.erate <- data.frame(bs.erate)

  ## Assign names to bs.erate
  for (j in 1:d.f) {
    names(bs.erate)[j] <- paste("E.rate", j, sep = "")
    names(bs.p.dist)[j] <- paste("p.dist", j, sep = "")
  }
  ## Rescale p.dist to be in [0,1]
  erate <- (erate - min(erate))/(max(erate)-min(erate))


  ### Transform the trait variables (t.data) to fill the gaps in the data for smooth plotting using seq
  ff <- c()
  tt <- c()
  if(ncom > n){
    myn <- n
  }else{
    myn <- ncom
  }

  t.data.trans <- matrix(NA, nrow = myn, ncol = ncol(t.data))
  for (j in 1:ncol(t.data)) {
    t.data.trans[,j] <- seq(from = range(t.data[,j])[1], to = range(t.data[,j])[2], length.out = myn)
  }
  t.data.trans <- `names<-`(data.frame(t.data.trans), names(t.data))

  #### Compute the B-splines of the transformed trait variables (t.data.trans) to get bs.traits.trans
  bs.traits.trans <- matrix(NA, nrow(t.data.trans), (ncol(t.data.trans) * d.f))
  for (j in 1:ncol(t.data.trans)) {
    for (i in 1:d.f) {
      bs.traits.trans[,(j - 1) * d.f + i] <- splines2::bSpline(t.data.trans[,j], d.f=d.f, degree=degree, intercept = TRUE)[,i]
    }
  }
  bs.traits.trans <- data.frame(bs.traits.trans)

  ## Assign names to bs.traits.trans
  for (i in 1:(ncol(bs.traits.trans)/d.f)) {
    for (j in 1:d.f) {
      names(bs.traits.trans)[(i - 1) * d.f + j] <- paste(names(t.data.trans)[i], j, sep = "")
    }
  }

  ### Transform the p.dist & erate variables to fill the gaps in the data for smooth plotting using seq
  p.dist.trans <- seq(from = range(p.dist)[1], to = range(p.dist)[2], length.out = myn)
  erate.trans <- seq(from = range(erate)[1], to = range(erate)[2], length.out = myn)

  ### Compute the B-splines of the transformed p.dist & erate variables (p.dist.trans & erate.trans) to get bs.p.dist.trans & bs.erate.trans
  bs.erate.trans <- matrix(NA, nrow=myn, ncol=d.f)
  bs.p.dist.trans <- matrix(NA, nrow=myn, ncol=d.f)
  for (i in 1:d.f) {
    bs.erate.trans[,i] <- splines2::bSpline(erate.trans, d.f=d.f, degree=degree, intercept = TRUE)[,i]
    bs.p.dist.trans[,i] <- splines2::bSpline(p.dist.trans, d.f=d.f, degree=degree, intercept = TRUE)[,i]
  }
  bs.erate.trans <- data.frame(bs.erate.trans)
  bs.p.dist.trans <- data.frame(bs.p.dist.trans)

  ## Assign names to bs.erate.trans and bs.p.dist.trans
  for (j in 1:d.f) {
    names(bs.erate.trans)[j] <- paste("E.rate", j, sep = "")
    names(bs.p.dist.trans)[j] <- paste("p.dist", j, sep = "")
  }

  ## Transformed values of all variables (trans.variables) and their B-splines (bs.variables.trans)
  bs.variables.trans <- cbind(bs.traits.trans, bs.p.dist.trans, bs.erate.trans)
  trans.variables <- cbind(t.data.trans, p.dist.trans, erate.trans)
  names(trans.variables) <- c(names(t.data.trans), "P.dist", "E.rate")


  #### Plot to see if the correct B-spline plots are output
  if(b.plots==TRUE){
    cols <- rep(c("red","blue","black","green"), 8*d.f)
    if(bsplines=="single"){
      grDevices::pdf(file = paste0(system.file("ms", package = "msco"), "/B-splines.curves_single.predictor.pdf"), height = 5, width = 5)
      plot(x=trans.variables[,1], y=bs.variables.trans[,(1+((1-1)*d.f))], type = "l", lty=1, lwd=2, col=cols[1], ylim=c(0,max(bs.variables.trans[,(1+((1-1)*d.f))])),
           xlab = "Trait value", ylab = "B-splines", main = "Trait variable")

      for(i in 2:4){
        graphics::lines(trans.variables[,1], bs.variables.trans[,i], col=cols[i], lwd=2, lty = i)
      }
      for (i in 1:4) {
        graphics::points(t.data[,1], bs.traits[,i], col=cols[i], pch=match(cols[i], cols))
      }
      graphics::text(0.1, 0.95, expression(paste(B["0,4"])), font=2)
      graphics::text(0.35, 0.52, expression(paste(B["1,4"])), font=2)
      graphics::text(0.7, 0.52, expression(paste(B["2,4"])), font=2)
      graphics::text(0.9, 0.95, expression(paste(B["3,4"])), font=2)

      grDevices::dev.off()
      print(noquote("Check msco's 'ms' folder in your R version's directory for a 'B-splines.curves_single.predictor.pdf' file."))
    }else{
      grDevices::pdf(file = paste0(system.file("ms", package = "msco"), "/B-splines.curves_all.predictors.pdf"), height = 8.27, width = 6)
      graphics::par(mar=c(4,4,2,0.5)+.1)
      graphics::par(mfcol=c(ceiling(((ncol(t.data)+2)/2)), 2))

      for (v.index in 1:ncol(trans.variables)) {
        plot(x=trans.variables[,v.index], y=bs.variables.trans[,(1+((v.index-1)*d.f))], type = "l", lty=1, lwd=2, col=cols[1], ylim=c(0,max(bs.variables.trans[,(1+((v.index-1)*d.f))])),
             xlab = "Trait value", ylab = "B-splines", main = paste(names(trans.variables)[v.index]))

        for(i in (((v.index-1)*d.f)+2):(((v.index-1)*d.f)+d.f)){
          graphics::lines(trans.variables[,v.index], bs.variables.trans[,i], col=cols[i], lwd=2, lty = i)
        }
        if(v.index < ncol(t.data) | v.index==ncol(t.data)){
          for (i in (1+((v.index-1)*d.f)):(((v.index-1)*d.f)+d.f)) {
            graphics::points(t.data[,v.index], bs.traits[,i], col=cols[i], pch=match(cols[i], cols))
          }
        }else if(v.index==(ncol(t.data)+1)){
          for (i in 1:d.f) {
            graphics::points(p.dist, bs.p.dist[,i], col=cols[i], pch=match(cols[i], cols))
          }
        }else if(v.index==(ncol(t.data)+2)){
          for (i in 1:d.f) {
            graphics::points(erate, bs.erate[,i], col=cols[i], pch=match(cols[i], cols))
          }
        }
        graphics::text(0.13, 0.9, expression(paste(B["0,4"])), font=2)
        graphics::text(0.35, 0.55, expression(paste(B["1,4"])), font=2)
        graphics::text(0.7, 0.55, expression(paste(B["2,4"])), font=2)
        graphics::text(0.87, 0.9, expression(paste(B["3,4"])), font=2)
      }
      grDevices::dev.off()
      print(noquote("Check msco's 'ms' folder in your R version's directory for a 'B-splines.curves_all.predictors.pdf' file."))
    }
  }

  ## Confirm if the plot of the sum of the B-splines is constant at y=1
  # plot(t.data[,1], rowSums(bs.traits[, seq(1,ncol(bs.traits),d.f)[1]:seq(d.f,ncol(bs.traits),d.f)[1]]), col=cols[1])
  # for (i in 2:((ncol(bs.traits))/d.f)) {
  #   points(t.data[,i], rowSums(bs.traits[, seq(1,ncol(bs.traits),d.f)[i]:seq(d.f,ncol(bs.traits),d.f)[i]]), col=cols[i])
  # }

  ######################################################################################################################
  ############################### END (for basis functions plots) ######################################################


  ######################################################################################################################
  ############################### B-spline regression (start) ##########################################################
  ######################################################################################################################

  ## Combine bs.traits.diff, bs.p.dist (b-splines of p.dist) and bs.erate (b-splines of erate)
  bs.variables.diff <- cbind(bs.traits.diff, bs.p.dist, bs.erate)

  ## Perform regression between jo and bs.variables.diff to get the regression coefficients
  bs.variables.diff <- data.frame(bs.variables.diff)
  if((metric %in% c("raw_prop", "Simpson_eqn", "Sorensen_eqn"))==TRUE){
    model <- suppressWarnings(glm2::glm2(jo ~ ., family=stats::quasibinomial(link="log"), data = bs.variables.diff,
                                         start = seq(start.range[1], start.range[2], length.out=(ncol(bs.variables.diff))+1)))
  }else if(metric=="raw" & (metric %in% c("raw_prop", "Simpson_eqn", "Sorensen_eqn"))!=TRUE & gbsm.model=="quasipoisson"){
    model <- suppressWarnings(glm2::glm2(jo ~ ., family=stats::quasipoisson(link="log"), data = bs.variables.diff,
                                         start = seq(start.range[1], start.range[2], length.out=(ncol(bs.variables.diff))+1)))
  }else if(metric=="raw" & (metric %in% c("raw_prop", "Simpson_eqn", "Sorensen_eqn"))!=TRUE & gbsm.model=="nb"){
    model <- suppressWarnings(MASS::glm.nb(jo ~ ., link = log, data = bs.variables.diff,
                                           start = seq(start.range[1], start.range[2], length.out=(ncol(bs.variables.diff))+1)))
  }else if(metric=="raw" & (metric %in% c("raw_prop", "Simpson_eqn", "Sorensen_eqn"))!=TRUE & (gbsm.model %in% c("quasipoisson", "nb"))!=TRUE){
    stop("Wrong 'gbsm.model' used for 'raw' version of joint occupancy. It must either be 'quasipoisson' or 'nb'.")
  }

  mysum <- summary(model)

  ##### Variance explained
  # pred.jo <- stats::predict(model, newdata = bs.variables.diff, type = "response", se = TRUE)
  pred.jo <- suppressWarnings(stats::predict.glm(model, newdata = bs.variables.diff, type = "response"))
  pred.jo <- as.numeric(pred.jo)
  var.expd2 <- (stats::cor(jo, pred.jo))^2

  if(scat.plot==TRUE){
    plot(jo, pred.jo, xlab="Joint occupancy", ylab="Predicted J. occ")
  }


  ## Obtain the regression coefficients excluding the intercept

  coeff <- model$coefficients[-1]
  coeff[which(is.na(coeff)==TRUE)] <- 0
  intercept <- model$coefficients[1]

  #####################################################################
  ############# Weighted sum of responses  ############################
  #####################################################################

  ########## traits (t.data)
  J_preds.traits <- matrix(NA, nrow=nrow(bs.traits), ncol=ncol(bs.traits))
  coeff.traits <- coeff[1:dim(bs.traits)[2]]
  for (i in 1:ncol(bs.traits)) {
    J_preds.traits[,i] <- coeff.traits[i]*bs.traits[,i] + intercept/((dim(t.data)[2]+2)*d.f)
  }
  J_preds.traits <- `names<-`(data.frame(J_preds.traits), names(bs.traits))

  # Sum "d.f" columns per variable (from J_preds.traits) to get "J_preds.traits.fin"
  names(J_preds.traits) <- gsub("[[:digit:]]", "", names(J_preds.traits))
  J_preds.traits.fin <- t(rowsum(t(J_preds.traits), group = colnames(J_preds.traits), na.rm = T, reorder=FALSE))

  ########## p.dist
  J_preds.p.d <- matrix(NA, nrow=nrow(bs.p.dist), ncol=ncol(bs.p.dist))
  coeff.p.d <- coeff[(dim(bs.traits)[2]+1):(dim(bs.traits)[2]+d.f)]
  for (i in 1:ncol(bs.p.dist)) {
    J_preds.p.d[,i] <- coeff.p.d[i]*bs.p.dist[,i] + intercept/((dim(t.data)[2]+2)*d.f)
  }
  J_preds.p.d <- `names<-`(data.frame(J_preds.p.d), names(bs.p.dist))

  # Sum "d.f" columns per variable (from J_preds.p.d) to get "J_preds.p.d.fin"
  names(J_preds.p.d) <- gsub("[[:digit:]]", "", names(J_preds.p.d))
  J_preds.p.d.fin <- t(rowsum(t(J_preds.p.d), group = colnames(J_preds.p.d), na.rm = TRUE, reorder=FALSE))

  ########## E.rate
  J_preds.er <- matrix(NA, nrow=nrow(bs.erate), ncol=ncol(bs.erate))
  coeff.er <- coeff[(dim(bs.traits)[2]+d.f+1):(dim(bs.variables.diff)[2])]
  for (i in 1:ncol(bs.erate)) {
    J_preds.er[,i] <- coeff.er[i]*bs.erate[,i] + intercept/((dim(t.data)[2]+2)*d.f)
  }
  J_preds.er <- `names<-`(data.frame(J_preds.er), names(bs.erate))

  # Sum "d.f" columns per variable (from J_preds.er) to get "J_preds.er.fin"
  names(J_preds.er) <- gsub("[[:digit:]]", "", names(J_preds.er))
  J_preds.er.fin <- t(rowsum(t(J_preds.er), group = colnames(J_preds.er), na.rm = TRUE, reorder=FALSE))


  #############################################################################
  ############# Weighted sum of seq-transformed responses #####################
  #############################################################################

  ##Product of coefficients with variables
  coeff.variables <- c(coeff.traits, coeff.p.d, coeff.er)
  J_preds.trans <- matrix(NA, nrow=nrow(bs.variables.trans), ncol=ncol(bs.variables.trans))
  for (i in 1:ncol(bs.variables.trans)) {
    J_preds.trans[,i] <- coeff.variables[i]*bs.variables.trans[,i] + intercept/((dim(t.data)[2]+2)*d.f)
  }
  J_preds.trans <- `names<-`(data.frame(J_preds.trans), names(bs.variables.trans))
  ## Sum "d.f" columns per variable (from J_preds.trans) to get "J_preds.trans.fin"
  names(J_preds.trans) <- gsub("[[:digit:]]", "", names(J_preds.trans))
  J_preds.trans.fin <- t(rowsum(t(J_preds.trans), group = colnames(J_preds.trans), na.rm = TRUE, reorder=FALSE))


  ####### Plots

  if(response.curves==TRUE){
    cols <- c("red","blue","black","green", "orange", "brown", "purple")
    limits <- range(J_preds.trans.fin)
    llim <- limits[1]
    ulim <- limits[2]
    if(ylabel==TRUE){
      if(leg==0){
        plot(x=trans.variables[,1], y=J_preds.trans.fin[,1], type = "l", lty=1, lwd=1.5, col=cols[1], ylim=c(llim, ulim),
             xlab = "Predictor value", ylab = "Wtd pred. value", main = noquote(paste("Order",order.jo)), cex=0.8, pch=1)
        for(i in 2:ncol(trans.variables)){
          graphics::lines(trans.variables[,i], J_preds.trans.fin[,i], col=cols[i], lwd=1.5, lty=i, cex=0.8, pch=i)
        }
        for(i in 1:ncol(t.data)){
          graphics::points(t.data[,i], J_preds.traits.fin[,i], col=cols[i], lwd=1, lty=i, cex=0.7, pch=i)
        }
        graphics::points(p.dist, J_preds.p.d.fin, col=cols[(ncol(t.data)+1)], lwd=1, pch=10, cex=0.6)
        graphics::points(erate, J_preds.er.fin, col=cols[(ncol(t.data)+2)], lwd=1, pch=18, cex=0.7)
      }else if(leg==1){
        plot(x=trans.variables[,1], y=J_preds.trans.fin[,1], type = "l", lty=1, lwd=1.5, col=cols[1], ylim=c(llim, (ulim + 0.2)),
             xlab = "Predictor value", ylab = "Wtd pred. value", main = noquote(paste("Order",order.jo)), cex=0.8, pch=1)
        for(i in 2:ncol(trans.variables)){
          graphics::lines(trans.variables[,i], J_preds.trans.fin[,i], col=cols[i], lwd=1.5, lty=i, cex=0.8, pch=i)
        }
        for(i in 1:ncol(t.data)){
          graphics::points(t.data[,i], J_preds.traits.fin[,i], col=cols[i], lwd=1, lty=i, cex=0.7, pch=i)
        }
        graphics::points(p.dist, J_preds.p.d.fin, col=cols[(ncol(t.data)+1)], lwd=1, pch=10, cex=0.6)
        graphics::points(erate, J_preds.er.fin, col=cols[(ncol(t.data)+2)], lwd=1, pch=18, cex=0.7)
        graphics::legend("top", legend = names(trans.variables), col = cols, lty=1:ncol(trans.variables), lwd=1.5,
                         pch = 1:(length(t.data)+2), bty = "n", cex = 0.8, ncol = 2)

      }else if(leg==2){
        plot(x=trans.variables[,1], y=J_preds.trans.fin[,1], type = "l", lty=1, lwd=1.5, col=cols[1], ylim=c(llim, (ulim + 2)),
             xlab = "Predictor value", ylab = "Wtd pred. value", main = noquote(paste("Order",order.jo)), cex=0.8, pch=1)
        for(i in 2:ncol(trans.variables)){
          graphics::lines(trans.variables[,i], J_preds.trans.fin[,i], col=cols[i], lwd=1.5, lty=i, cex=0.8, pch=i)
        }
        for(i in 1:ncol(t.data)){
          graphics::points(t.data[,i], J_preds.traits.fin[,i], col=cols[i], lwd=1, lty=i, cex=0.7, pch=i)
        }
        graphics::points(p.dist, J_preds.p.d.fin, col=cols[(ncol(t.data)+1)], lwd=1, pch=10, cex=0.6)
        graphics::points(erate, J_preds.er.fin, col=cols[(ncol(t.data)+2)], lwd=1, pch=18, cex=0.7)
        graphics::legend("top", legend = names(trans.variables), col = cols, lty=1:ncol(trans.variables), lwd=1.5,
                         pch = 1:(length(t.data)+2), bty = "n", cex = 0.8, ncol = 2)
      }
    }else{
      if(leg==0){
        plot(x=trans.variables[,1], y=J_preds.trans.fin[,1], type = "l", lty=1, lwd=1.5, col=cols[1], ylim=c(llim, ulim),
             xlab = "Predictor value", ylab = " ", main = noquote(paste("Order",order.jo)), cex=0.8, pch=1)
        for(i in 2:ncol(trans.variables)){
          graphics::lines(trans.variables[,i], J_preds.trans.fin[,i], col=cols[i], lwd=1.5, lty=i, cex=0.8, pch=i)
        }
        for(i in 1:ncol(t.data)){
          graphics::points(t.data[,i], J_preds.traits.fin[,i], col=cols[i], lwd=1, lty=i, cex=0.7, pch=i)
        }
        graphics::points(p.dist, J_preds.p.d.fin, col=cols[(ncol(t.data)+1)], lwd=1, pch=10, cex=0.6)
        graphics::points(erate, J_preds.er.fin, col=cols[(ncol(t.data)+2)], lwd=1, pch=18, cex=0.7)
      }else if(leg==1){
        plot(x=trans.variables[,1], y=J_preds.trans.fin[,1], type = "l", lty=1, lwd=1.5, col=cols[1], ylim=c(llim, (ulim + 0.2)),
             xlab = "Predictor value", ylab = " ", main = noquote(paste("Order",order.jo)), cex=0.8, pch=1)
        for(i in 2:ncol(trans.variables)){
          graphics::lines(trans.variables[,i], J_preds.trans.fin[,i], col=cols[i], lwd=1.5, lty=i, cex=0.8, pch=i)
        }
        for(i in 1:ncol(t.data)){
          graphics::points(t.data[,i], J_preds.traits.fin[,i], col=cols[i], lwd=1, lty=i, cex=0.7, pch=i)
        }
        graphics::points(p.dist, J_preds.p.d.fin, col=cols[(ncol(t.data)+1)], lwd=1, pch=10, cex=0.6)
        graphics::points(erate, J_preds.er.fin, col=cols[(ncol(t.data)+2)], lwd=1, pch=18, cex=0.7)
        graphics::legend("top", legend = names(trans.variables), col = cols, lty=1:ncol(trans.variables), lwd=1.5,
                         pch = 1:(length(t.data)+2), bty = "n", cex = 0.8, ncol = 2)

      }else if(leg==2){
        plot(x=trans.variables[,1], y=J_preds.trans.fin[,1], type = "l", lty=1, lwd=1.5, col=cols[1], ylim=c(llim, (ulim + 2)),
             xlab = "Predictor value", ylab = " ", main = noquote(paste("Order",order.jo)), cex=0.8, pch=1)
        for(i in 2:ncol(trans.variables)){
          graphics::lines(trans.variables[,i], J_preds.trans.fin[,i], col=cols[i], lwd=1.5, lty=i, cex=0.8, pch=i)
        }
        for(i in 1:ncol(t.data)){
          graphics::points(t.data[,i], J_preds.traits.fin[,i], col=cols[i], lwd=1, lty=i, cex=0.7, pch=i)
        }
        graphics::points(p.dist, J_preds.p.d.fin, col=cols[(ncol(t.data)+1)], lwd=1, pch=10, cex=0.6)
        graphics::points(erate, J_preds.er.fin, col=cols[(ncol(t.data)+2)], lwd=1, pch=18, cex=0.7)
        graphics::legend("top", legend = names(trans.variables), col = cols, lty=1:ncol(trans.variables), lwd=1.5,
                         pch = 1:(length(t.data)+2), bty = "n", cex = 0.8, ncol = 2)
      }
    }
  }

  gbs <- list()
  gbs$Predictors <- trans.variables
  gbs$Responses <- as.data.frame(J_preds.trans.fin)
  gbs$order.jo <- order.jo
  gbs$coeff <- coeff
  gbs$glm_obj <- model
  gbs$j.occs <- jo
  gbs$pred.j.occs <- pred.jo
  gbs$bs_pred <- bs.variables.diff
  gbs$start.range <- start.range
  gbs$var.expld <- var.expd2
  gbs$summary <- mysum
  class(gbs) <- "gbsm"
  return(gbs)
}

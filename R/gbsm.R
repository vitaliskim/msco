#' A generalised B-spline modelling for a set of neutral and trait-based variables
#'
#' This function implements the generalised B-spline model (*sensu* Lagat *et al.,* 2021b)
#'  for dissecting the effects of neutral encounter versus functional traits on multi-order
#'   species interactions and co-occurrence. Generalized linear model
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
#'  `Simpson_eqn` for Simpson equivalent, `Sorensen_eqn` for Sorensen equivalent, and `raw` for the
#'   raw form of index without rescaling.
#' @param n Number of samples for which the joint occupancy is computed. These samples are non-overlapping.
#'  I.e., sampling is done without replacement. If the total number of combinations of `i` species chosen
#'   from the total species pool `m`, i.e. `choose(m,i)`, is less than this value (`n`), `choose(m,i)` is
#'    used as the (maximum) number of samples one can set. Otherwise sampling without replacement is
#'     performed to select just the `n` samples.
#' @param b.plots Boolean value indicating if B-spline basis functions should be plotted.
#' @param bsplines This parameter indicates if a single or all B-spline curves should be plotted.
#'  If `b.plots=TRUE` and `bsplines="single"`, the B-spline curves for the first predictor in
#'   `t.data` will be plotted. Any other value for `bsplines` (other than `"single"`) results in the
#'    B-spline curves for all predictors being plotted.
#' @param response.curves A boolean value indicating if all response curves should be plotted.
#' @param leg Boolean value indicating if the legend of the gbsm outputs should be included in the plots. This
#'  parameter is added to help control the appearance of plots in \link[msco]{gbsm_m.orders} function.
#' @param scat.plot Boolean value indicating if scatter plots between joint occupancy and its predicted
#'  values should be plotted.
#' @param start Starting values for glm regression.
#'
#' @return `gbsm` function returns a list containing the following outputs:
#' \item{`order.jo`}{Order of joint occupancy}
#' \item{`Predictors`}{Predictor variables used in GLM regression with binomial variance
#'  distribution function and log link function.}
#' \item{`Responses`}{Response variables from GLM regression with binomial variance distribution
#'  function and log link function.}
#' \item{`coeff`}{Coefficients of the generalized linear model used.}
#' \item{`glm_obj`}{Generalized linear model used.}
#' \item{`j.occs`}{Observed joint occupancies.}
#' \item{`pred.j.occs`}{Predicted joint occupancies.}
#' \item{`bs_pred`}{B-spline-transformed `Predictors`.}
#' \item{`start`}{Starting values for the generalized linear model used.}
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
#'  \item{Lagat, V. K., Latombe, G. and Hui, C. (2021b). *Dissecting the effects of
#'   neutral encounter versus functional traits on multi-order species interactions
#'    and co-occurrence with generalised B-spline modelling*. DOI: `<To be added>`.}
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
#'  my.gbsm <- msco::gbsm(s.data, t.data, p.d.mat, metric = "Simpson_eqn",
#'   d.f=4, order.jo=3, degree=3, n=1000, b.plots=TRUE, scat.plot=TRUE,
#'    bsplines="single", response.curves=TRUE, leg=1,
#'     start=seq(-0.1, 0, length.out=(ncol(t.data)+2)*4+1))
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
                 bsplines="single", scat.plot=TRUE, response.curves=TRUE, leg=1, start=seq(-0.1, 0, length.out=(ncol(t.data)+2)*4+1)){

  if(class(t.data)!="data.frame"){
    t.data <- as.data.frame(t.data)
  }

  # if (length(setdiff(sapply(t.data, class), c("factor", "numeric"))) > 0) {
  #   stop("Trait variables must be factor or numeric")
  # }
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

  #### Compute the B-splines of the original trait variables (t.data) to get bt.data
  bt.data <- matrix(NA, nrow(t.data), (ncol(t.data) * d.f))
  for (j in 1:ncol(t.data)) {
    for (i in 1:d.f) {
      bt.data[,(j - 1) * d.f + i] <- splines2::bSpline(t.data[,j], d.f=d.f, degree=degree, intercept = TRUE)[,i]
    }
  }
  bt.data <- data.frame(bt.data)

  ## Assign names to bt.data
  for (i in 1:(ncol(bt.data)/d.f)) {
    for (j in 1:d.f) {
      names(bt.data)[(i - 1) * d.f + j] <- paste(names(t.data)[i], j, sep = "")
    }
  }
  bt.data <- `rownames<-`(bt.data, rownames(t.data))

  ## Compute the differences on transformed values (bt.data) using SD, to get t.mat. Simultaneously compute p.dist and encounter rate
  order <- 1:nrow(s.data) ## Possible joint occupancy orders
  sn <- dim(s.data)[1] ## Total number of species
  ncom <- choose(sn,order.jo) ## Total number of the combinations of "order.jo" species chosen from sn

  if(ncom > n){ ## Use MCMC to sample n species combinations if the total combinations > the chosen sample size, n.
    jo <- rep(NA, n)
    t.mat <- matrix(NA, nrow = n, ncol = ncol(bt.data))
    erate <- rep(NA, n)
    p.dist <- rep(NA, n)
    for (j in 1:n) {
      sam <- sample(1:sn, order.jo, replace = FALSE)
      t.mat[j,] <- apply(bt.data[sam,], 2, stats::sd)

      #### jo values for chosen combination of species
      if(metric=="raw"){
        jo[j] <- msco::j.occ(s.data[sam,], order = order.jo)$jo.val
      }else if(metric=="Simpson_eqn"){
        jo[j] <- (msco::j.occ(s.data[sam,], order = order.jo)$jo.val)/min(rowSums(s.data[sam,]))
      }else if(metric=="Sorensen_eqn"){
        jo[j] <- (msco::j.occ(s.data[sam,], order = order.jo)$jo.val)/mean(rowSums(s.data[sam,]))
      }
      ## p.dist
      p.dist[j] <- mean(p.d.mat[t(utils::combn(sort(sam), 2))])

      ## Encounter rate
      nsam <- length(sam)
      er <- rep(NA, nsam)
      for (i in 1:nsam) {
        er[i] <- rowSums(s.data[sam[i],])
      }
      erate[j] <- (prod(er))/((ncol(s.data))^nsam)
    }
    if((metric %in% c("raw", "Simpson_eqn", "Sorensen_eqn"))!=TRUE){
      stop("Wrong option for the joint occupancy metric provided. It must either be 'raw', 'Simpson_eqn', or 'Sorensen_eqn'.")
    }
  }else{ ## Otherwise use the combinations as they are without sampling n of them.
    jo <- rep(NA, ncom)
    t.mat <- matrix(NA, nrow = ncom, ncol = ncol(bt.data))
    erate <- rep(NA, ncom)
    p.dist <- rep(NA, ncom)
    com <- utils::combn(sn, order.jo)
    for (j in 1:ncom) {
      t.mat[j,] <- apply(bt.data[com[,j],], 2, stats::sd)

      #### jo values for chosen combination of species
      if(metric=="raw"){
        jo[j] <- msco::j.occ(s.data[com[,j],], order = order.jo)$jo.val
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
      erate[j] <- (prod(er))/((ncol(s.data))^nsam)
    }
    if((metric %in% c("raw", "Simpson_eqn", "Sorensen_eqn"))!=TRUE){
      stop("Wrong option for the joint occupancy metric provided. It must either be 'raw', 'Simpson_eqn', or 'Sorensen_eqn'.")
    }
  }

  p.dist <- (p.dist-min(p.dist))/(max(p.dist)-min(p.dist))
  erate <- (erate-min(erate))/(max(erate)-min(erate))

  ## B-splines of p.dist & erate
  bt.erate <- matrix(NA, nrow=length(erate), ncol=d.f)
  bt.p.dist <- matrix(NA, nrow=length(p.dist), ncol=d.f)
  for (i in 1:d.f) {
    bt.erate[,i] <- splines2::bSpline(erate, d.f=d.f, degree=degree, intercept = TRUE)[,i]
    bt.p.dist[,i] <- splines2::bSpline(p.dist, d.f=d.f, degree=degree, intercept = TRUE)[,i]
  }
  bt.erate <- data.frame(bt.erate)
  bt.p.dist <- data.frame(bt.p.dist)

  ## Assign names to bt.erate & bt.p.dist
  for (j in 1:d.f) {
    names(bt.erate)[j] <- paste("E.rate", j, sep = "")
    names(bt.p.dist)[j] <- paste("p.dist", j, sep = "")
  }

  ### Transform the trait variables (t.data) to fill the gaps in the data for smooth plotting
  ff <- c()
  tt <- c()
  if(ncom > n){
    myn <- n
  }else{
    myn <- ncom
  }

  t.trans <- matrix(NA, nrow = myn, ncol = ncol(t.data))
  for (j in 1:ncol(t.data)) {
    t.trans[,j] <- seq(from = range(t.data[,j])[1], to = range(t.data[,j])[2], length.out = myn)
  }
  t.trans <- `names<-`(data.frame(t.trans), names(t.data))

  #### Compute the B-splines of the transformed trait variables (t.trans) to get bt.trans
  bt.trans <- matrix(NA, nrow(t.trans), (ncol(t.trans) * d.f))
  for (j in 1:ncol(t.trans)) {
    for (i in 1:d.f) {
      bt.trans[,(j - 1) * d.f + i] <- splines2::bSpline(t.trans[,j], d.f=d.f, degree=degree, intercept = TRUE)[,i]
    }
  }
  bt.trans <- data.frame(bt.trans)

  ## Assign names to bt.trans
  for (i in 1:(ncol(bt.trans)/d.f)) {
    for (j in 1:d.f) {
      names(bt.trans)[(i - 1) * d.f + j] <- paste(names(t.trans)[i], j, sep = "")
    }
  }

  ### Transform the p.dist & erate variables to fill the gaps in the data for smooth plotting
  pd.trans <- seq(from = range(p.dist)[1], to = range(p.dist)[2], length.out = myn)
  e.trans <- seq(from = range(erate)[1], to = range(erate)[2], length.out = myn)

  ### Compute the B-splines of the transformed p.dist & erate variables (pd.trans & e.trans) to get bt.pdtrans & bt.etrans
  bt.etrans <- matrix(NA, nrow=myn, ncol=d.f)
  bt.pdtrans <- matrix(NA, nrow=myn, ncol=d.f)
  for (i in 1:d.f) {
    bt.etrans[,i] <- splines2::bSpline(e.trans, d.f=d.f, degree=degree, intercept = TRUE)[,i]
    bt.pdtrans[,i] <- splines2::bSpline(pd.trans, d.f=d.f, degree=degree, intercept = TRUE)[,i]
  }
  bt.etrans <- data.frame(bt.etrans)
  bt.pdtrans <- data.frame(bt.pdtrans)

  ## Assign names to bt.etrans and bt.pdtrans
  for (j in 1:d.f) {
    names(bt.etrans)[j] <- paste("E.rate", j, sep = "")
    names(bt.pdtrans)[j] <- paste("p.dist", j, sep = "")
  }

  ## Transformed values of all variables (t.trans) and their B-splines (bt.trans)
  bt.trans <- cbind(bt.trans, bt.pdtrans, bt.etrans)
  t.transs <- cbind(t.trans, pd.trans, e.trans)
  names(t.transs) <- c(names(t.trans), "P.dist", "E.rate")
  t.trans <- t.transs

  #### Plot to see if the correct B-spline plots are output
  if(b.plots==TRUE){
    cols <- rep(c("red","blue","black","green"), 8*d.f)
    if(bsplines=="single"){
      grDevices::pdf(file = paste0(system.file("ms", package = "msco"), "/B-splines.curves_single.predictor.pdf"), height = 5, width = 5)
      plot(x=t.trans[,1], y=bt.trans[,(1+((1-1)*d.f))], type = "l", lty=1, lwd=2, col=cols[1], ylim=c(0,max(bt.trans[,(1+((1-1)*d.f))])),
           xlab = "Trait variable", ylab = "B-spline curves", main = paste(names(t.trans)[1]))

      for(i in 2:4){
        graphics::lines(t.trans[,1], bt.trans[,i], col=cols[i], lwd=2, lty = i)
      }
      for (i in 1:4) {
        graphics::points(t.data[,1], bt.data[,i], col=cols[i], pch=match(cols[i], cols))
      }
      graphics::text(0.1, 0.95, expression(paste(B["0,4"])), font=2)
      graphics::text(0.35, 0.52, expression(paste(B["1,4"])), font=2)
      graphics::text(0.7, 0.52, expression(paste(B["2,4"])), font=2)
      graphics::text(0.9, 0.95, expression(paste(B["3,4"])), font=2)

      grDevices::dev.off()
      base::system(paste0('open "', paste0(system.file("ms", package = "msco"), "/B-splines.curves_single.predictor.pdf"), '"'))
    }else{
      grDevices::pdf(file = paste0(system.file("ms", package = "msco"), "/B-splines.curves_all.predictors.pdf"), height = 8.27, width = 6)
      graphics::par(mar=c(4,4,2,0.5)+.1)
      graphics::par(mfcol=c(ceiling(((ncol(t.data)+2)/2)), 2))

      for (v.index in 1:ncol(t.trans)) {
        plot(x=t.trans[,v.index], y=bt.trans[,(1+((v.index-1)*d.f))], type = "l", lty=1, lwd=2, col=cols[1], ylim=c(0,max(bt.trans[,(1+((v.index-1)*d.f))])),
             xlab = "Trait variable", ylab = "B-spline curves", main = paste(names(t.trans)[v.index]))

        for(i in (((v.index-1)*d.f)+2):(((v.index-1)*d.f)+d.f)){
          graphics::lines(t.trans[,v.index], bt.trans[,i], col=cols[i], lwd=2, lty = i)
        }
        if(v.index < ncol(t.data) | v.index==ncol(t.data)){
          for (i in (1+((v.index-1)*d.f)):(((v.index-1)*d.f)+d.f)) {
            graphics::points(t.data[,v.index], bt.data[,i], col=cols[i], pch=match(cols[i], cols))
          }
        }else if(v.index==(ncol(t.data)+1)){
          for (i in 1:d.f) {
            graphics::points(p.dist, bt.p.dist[,i], col=cols[i], pch=match(cols[i], cols))
          }
        }else if(v.index==(ncol(t.data)+2)){
          for (i in 1:d.f) {
            graphics::points(erate, bt.erate[,i], col=cols[i], pch=match(cols[i], cols))
          }
        }
        graphics::text(0.13, 0.9, expression(paste(B["0,4"])), font=2)
        graphics::text(0.35, 0.55, expression(paste(B["1,4"])), font=2)
        graphics::text(0.7, 0.55, expression(paste(B["2,4"])), font=2)
        graphics::text(0.87, 0.9, expression(paste(B["3,4"])), font=2)
      }
      grDevices::dev.off()
      base::system(paste0('open "', paste0(system.file("ms", package = "msco"), "/B-splines.curves_all.predictors.pdf"), '"'))
    }
  }

  ## Confirm if the plot of the sum of the B-splines is constant at y=1
  # plot(t.data[,1], rowSums(bt.data[, seq(1,ncol(bt.data),d.f)[1]:seq(d.f,ncol(bt.data),d.f)[1]]), col=cols[1])
  # for (i in 2:((ncol(bt.data))/d.f)) {
  #   points(t.data[,i], rowSums(bt.data[, seq(1,ncol(bt.data),d.f)[i]:seq(d.f,ncol(bt.data),d.f)[i]]), col=cols[i])
  # }

  ######################################################################################################################
  ############################### END (for basis functions plots) ######################################################


  ######################################################################################################################
  ############################### B-spline regression (start) ##########################################################
  ######################################################################################################################

  ## Assign t.mat (trait differences of B-splines of t.data) column names
  t.mat <- `colnames<-`(t.mat, colnames(bt.data))
  t.mat <- data.frame(t.mat)

  ## Combine t.mat, p. distance and encounter rate
  t.mat <- cbind(t.mat, bt.p.dist, bt.erate)

  ## Rescale Xi values to be in [0,1] interval for raw jo
  if(metric=="raw"){
    jo <- (jo-min(jo))/(max(jo)-min(jo))
  }

  ## Perform regression between Xi and t.mat to get the regression coefficients
  t.mat <- data.frame(t.mat)
  model <- suppressWarnings(glm2::glm2(jo ~ ., family=stats::binomial(link="log"), data = t.mat, start = start))
  mysum <- summary(model)

  ##### Variance explained
  # pred.jo <- stats::predict(model, newdata = t.mat, type = "response", se = TRUE)
  pred.jo <- suppressWarnings(stats::predict.glm(model, newdata = t.mat, type = "response"))
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

  ########## t.data

  ## Product of the Coeffs with originally transformed matrix (bt.data)
  ott.mat <- matrix(NA, nrow=nrow(bt.data), ncol=ncol(bt.data))
  ct <- coeff[1:dim(bt.data)[2]]
  for (i in 1:ncol(bt.data)) {
    ott.mat[,i] <- ct[i]*bt.data[,i]
  }
  ott.mat <- `names<-`(data.frame(ott.mat), names(bt.data))
  ## Sum "d.f" columns per variable (from ott.mat) to get "ofdata"
  names(ott.mat) <- gsub("[[:digit:]]", "", names(ott.mat))
  ofdata <- t(rowsum(t(ott.mat), group = colnames(ott.mat), na.rm = T, reorder=FALSE))

  ########## p.dist
  ott.dis <- matrix(NA, nrow=nrow(bt.p.dist), ncol=ncol(bt.p.dist))
  cdis <- coeff[(dim(bt.data)[2]+1):(dim(bt.data)[2]+d.f)]

  for (i in 1:ncol(bt.p.dist)) {
    ott.dis[,i] <- cdis[i]*bt.p.dist[,i]
  }
  ott.dis <- `names<-`(data.frame(ott.dis), names(bt.p.dist))
  ## Sum "d.f" columns per variable (from ott.dis) to get "ofdist"
  names(ott.dis) <- gsub("[[:digit:]]", "", names(ott.dis))
  ofdist <- t(rowsum(t(ott.dis), group = colnames(ott.dis), na.rm = TRUE, reorder=FALSE))

  ########## E.rate
  ott.rat <- matrix(NA, nrow=nrow(bt.erate), ncol=ncol(bt.erate))
  ce <- coeff[(dim(bt.data)[2]+d.f+1):(dim(t.mat)[2])]
  for (i in 1:ncol(bt.erate)) {
    ott.rat[,i] <- ce[i]*bt.erate[,i]
  }
  ott.rat <- `names<-`(data.frame(ott.rat), names(bt.erate))
  ## Sum "d.f" columns per variable (from ott.mat) to get "ofdata"
  names(ott.rat) <- gsub("[[:digit:]]", "", names(ott.rat))
  ofrate <- t(rowsum(t(ott.rat), group = colnames(ott.rat), na.rm = TRUE, reorder=FALSE))


  #############################################################################
  ############# Weighted sum of transformed (using seq) responses #############
  #############################################################################


  coeff <- c(ct, cdis, ce)
  ott.trans <- matrix(NA, nrow=nrow(bt.trans), ncol=ncol(bt.trans))
  for (i in 1:ncol(bt.trans)) {
    ott.trans[,i] <- coeff[i]*bt.trans[,i]
  }
  ott.trans <- `names<-`(data.frame(ott.trans), names(bt.trans))
  ## Sum "d.f" columns per variable (from ott.trans) to get "oftrans"
  names(ott.trans) <- gsub("[[:digit:]]", "", names(ott.trans))
  oftrans <- t(rowsum(t(ott.trans), group = colnames(ott.trans), na.rm = TRUE, reorder=FALSE))

  ####### Plots

  if(response.curves==TRUE){
    cols <- c("red","blue","black","green", "orange", "brown", "purple")
    limits <- range(oftrans)
    llim <- limits[1]
    ulim <- limits[2]
    if(leg==0){
      plot(x=t.trans[,1], y=oftrans[,1], type = "l", lty=1, lwd=1.5, col=cols[1], ylim=c(llim, ulim),
           xlab = "Rescaled Range", ylab = "Effect on J. occ.", main = noquote(paste("Order",order.jo)), cex=0.8, pch=1)
      for(i in 2:ncol(t.trans)){
        graphics::lines(t.trans[,i], oftrans[,i], col=cols[i], lwd=1.5, lty=i, cex=0.8, pch=i)
      }
      for(i in 1:ncol(t.data)){
        graphics::points(t.data[,i], ofdata[,i], col=cols[i], lwd=1, lty=i, cex=0.7, pch=i)
      }
      graphics::points(p.dist, ofdist, col=cols[(ncol(t.data)+1)], lwd=1, pch=10, cex=0.6)
      graphics::points(erate, ofrate, col=cols[(ncol(t.data)+2)], lwd=1, pch=18, cex=0.7)
    }else if(leg==1){
      plot(x=t.trans[,1], y=oftrans[,1], type = "l", lty=1, lwd=1.5, col=cols[1], ylim=c(llim, (ulim + 0.2)),
           xlab = "Rescaled Range", ylab = "Effect on J. occ.", main = noquote(paste("Order",order.jo)), cex=0.8, pch=1)
      for(i in 2:ncol(t.trans)){
        graphics::lines(t.trans[,i], oftrans[,i], col=cols[i], lwd=1.5, lty=i, cex=0.8, pch=i)
      }
      for(i in 1:ncol(t.data)){
        graphics::points(t.data[,i], ofdata[,i], col=cols[i], lwd=1, lty=i, cex=0.7, pch=i)
      }
      graphics::points(p.dist, ofdist, col=cols[(ncol(t.data)+1)], lwd=1, pch=10, cex=0.6)
      graphics::points(erate, ofrate, col=cols[(ncol(t.data)+2)], lwd=1, pch=18, cex=0.7)
      graphics::legend("top", legend = names(t.trans), col = cols, lty=1:ncol(t.trans), lwd=1.5,
                       pch = 1:(length(t.data)+2), bty = "n", cex = 0.8, ncol = 2)

    }else if(leg==2){
      plot(x=t.trans[,1], y=oftrans[,1], type = "l", lty=1, lwd=1.5, col=cols[1], ylim=c(llim, (ulim + 2)),
           xlab = "Rescaled Range", ylab = "Effect on J. occ.", main = noquote(paste("Order",order.jo)), cex=0.8, pch=1)
      for(i in 2:ncol(t.trans)){
        graphics::lines(t.trans[,i], oftrans[,i], col=cols[i], lwd=1.5, lty=i, cex=0.8, pch=i)
      }
      for(i in 1:ncol(t.data)){
        graphics::points(t.data[,i], ofdata[,i], col=cols[i], lwd=1, lty=i, cex=0.7, pch=i)
      }
      graphics::points(p.dist, ofdist, col=cols[(ncol(t.data)+1)], lwd=1, pch=10, cex=0.6)
      graphics::points(erate, ofrate, col=cols[(ncol(t.data)+2)], lwd=1, pch=18, cex=0.7)
      graphics::legend("top", legend = names(t.trans), col = cols, lty=1:ncol(t.trans), lwd=1.5,
                       pch = 1:(length(t.data)+2), bty = "n", cex = 0.8, ncol = 2)
    }
  }

  gbs <- list()
  gbs$Predictors <- t.trans
  gbs$Responses <- as.data.frame(oftrans)
  gbs$order.jo <- order.jo
  gbs$coeff <- coeff
  gbs$glm_obj <- model
  gbs$j.occs <- jo
  gbs$pred.j.occs <- pred.jo
  gbs$bs_pred <- t.mat
  gbs$start <- start
  gbs$var.expld <- var.expd2
  gbs$summary <- mysum
  class(gbs) <- "gbsm"
  return(gbs)
}

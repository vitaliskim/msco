#' Predictor's contribution and model performance assessment from the results on multiple orders of joint occupancy
#'
#' This function implements the generalised B-spline model (gbsm; *sensu* Lagat et al., 2021b) for dissecting the
#'  effects of random encounter versus functional trait mismatching on multi-species co-occurrence and interference.
#'   Unlike \link[msco]{gbsm} that performs gbsm for a single order of species, \link[msco]{gbsm_m.orders} takes
#'    into account multiple orders of joint occupancy. In particular: for multiple joint occupancy orders, this
#'     function computes:
#' * each predictor's contribution to the explained variation in joint occupancy,
#' * the goodness-of-fit and model performance from cross-validation, and
#' * plots the:
#'      + response curves,
#'      + scatter plots (between the observed and predicted joint occupancy values),
#'      + histograms of the joint occupancy frequency distribution, and
#'      + model performance plots.
#'
#' @param orders Specific number of species for which the joint occupancy is computed.
#' @param degree As for \link[msco]{gbsm}.
#' @param n As for \link[msco]{gbsm}.
#' @param metric As for \link[msco]{gbsm}.
#' @param gbsm.model As for \link[msco]{gbsm}.
#' @param d.f As for \link[msco]{gbsm}.
#' @param k As for \link[msco]{cross_valid}.
#' @param p As for \link[msco]{cross_valid}.
#' @param type As for \link[msco]{cross_valid}.
#' @param s.data A species-by-site presence/absence `data.frame` with entries indicating
#' occurrence (1) and non-occurrence (0) of species in a site.
#' @param t.data A `data.frame` with traits as columns and species as rows. The species must be the same as
#'  in `s.data`.
#' @param p.d.mat A symmetric `matrix` with `dimnames` as species and entries indicating the
#'  phylogenetic distance between any two of them (species).
#' @param start As for \link[msco]{gbsm}.
#' @param scat.plots Boolean value indicating if scatter plots between joint occupancy and its predicted
#'  values should be plotted.
#' @param response.curves A boolean value indicating if all response curves for all joint occupancy
#'  orders (`jo.orders`) should be plotted.
#' @param j.occs.distrbn A boolean value indicating if the histograms of the frequency distribution of
#'  observed joint occupancy should be output.
#' @param mp.plots A boolean value indicating if the model performance plots should be output.
#' @return
#' `gbsm_m.orders` function returns a list containing the following outputs:
#'  * `jo.orders`: &nbsp; A set of joint occupancy orders.
#'  * `contrbn_table`: &nbsp; A `list` of `data.frame`s consisting of:
#'     + `predictor`: &nbsp; A column of predictors.
#'     + `var.expld_M1`: &nbsp; A column of goodness-of-fit (I.e., the Pearson's \eqn{r^2}
#'      between the observed and predicted values of joint occupancy when all predictors
#'       are used in the model.
#'     + `var.expld_M2`: &nbsp; The Pearson's \eqn{r^2} between the observed and the
#'      predicted values of joint occupancy when all predictors except the predictor
#'       whose contribution is to be determined, are used in the model.
#'     + `contribution`: &nbsp; Each predictor's proportion of contribution in
#'      explaining joint occupancy. This value is given by:
#'
#'          `contribution` = \eqn{\frac{var.expld_{M1} - var.expld_{M2}}{var.expld_{M1}}}
#'
#'  * `model.validation.table`: &nbsp; A `data.frame` with:
#'     + `orders`: &nbsp; Orders of joint occupancy used.
#'     + `Rsquared_gf`: &nbsp; Goodness-of-fit of the model. I.e., it is the
#'      Pearson's **\eqn{r^2}** between the observed and predicted values of
#'       joint occupancy, for different orders.
#'     + `Rsquared_cv`: &nbsp; Model performance from cross-validation.
#'
#'  * `metric`: &nbsp; As for \link[msco]{gbsm}.
#'  * `d.f`: &nbsp; As for \link[msco]{gbsm}.
#'  * `n`: &nbsp; As for \link[msco]{gbsm}.
#'  * `degree`: &nbsp; As for \link[msco]{gbsm}.
#'  * `jo.orders`: &nbsp; Orders of joint occupancy used.
#' @references
#' \enumerate{
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
#'   encounter versus functional trait mismatching on multi-species co-occurrence and
#'    interference with generalised B-spline modelling*. DOI: `<To be added>`.}
#' }
#' @examples
#' \dontrun{
#'  my.path <- system.file("extdata/gsmdat", package = "msco")
#'  setwd(my.path)
#'  s.data <- get(load("s.data.csv")) ## Species-by-site matrix
#'  t.data <- get(load("t.data.csv")) ## Species-by-Trait matrix
#'  p.d.mat <- get(load("p.d.mat.csv")) ## Species-by-species phylogenetic distance matrix
#'
#'
#'  RNGkind(sample.kind = "Rejection")
#'  set.seed(1)
#'  jp <- msco::gbsm_m.orders(s.data, t.data, p.d.mat, gbsm.model,
#'   metric="Simpson_eqn", orders = c(3:5, 8, 10, 15, 20), d.f=4,
#'    degree=3, n=1000, k=5, p=0.8, type="k-fold", scat.plots=TRUE,
#'     response.curves=TRUE, j.occs.distrbn=TRUE, mp.plots=TRUE,
#'      start=seq(-0.1, 0, length.out=(ncol(t.data)+2)*4+1))
#'
#'  jp$contbn_table[[1]]
#'  jp$model.validation.table
#'  jp$jo.orders
#'
#'  ## Close the open plots.gbsm.pdf file before running the 2nd example
#'  RNGkind(sample.kind = "Rejection")
#'  set.seed(1)
#'  jp2 <- msco::gbsm_m.orders(s.data, t.data, p.d.mat, gbsm.model,
#'   metric="Sorensen_eqn", orders = c(3:5, 8, 10, 15, 20), d.f=4,
#'    degree=3, n=1000, k=5, p=0.8, type="k-fold", scat.plots=TRUE,
#'     response.curves=TRUE, j.occs.distrbn=TRUE, mp.plots=TRUE,
#'      start=seq(-0.1, 0, length.out=(ncol(t.data)+2)*4+1))
#'
#'  jp2$contbn_table[[1]]
#'  jp2$model.validation.table
#'  jp2$jo.orders
#'
#' ## Close the open plots.gbsm.pdf file before running the 3rd example
#'  RNGkind(sample.kind = "Rejection")
#'  set.seed(1)
#'  jp3 <- msco::gbsm_m.orders(s.data, t.data, p.d.mat, gbsm.model,
#'   metric="raw_prop", orders = c(3:5, 8, 10, 15, 20), d.f=4,
#'    degree=3, n=1000, k=5, p=0.8, type="k-fold", scat.plots=TRUE,
#'     response.curves=TRUE, j.occs.distrbn=TRUE, mp.plots=TRUE,
#'      start=seq(-0.1, 0, length.out=(ncol(t.data)+2)*4+1))
#'
#'  jp3$contbn_table[[1]]
#'  jp3$model.validation.table
#'  jp3$jo.orders
#'
#' ## Close the open plots.gbsm.pdf file before running the 3rd example
#'  RNGkind(sample.kind = "Rejection")
#'  set.seed(1)
#'  jp4 <- msco::gbsm_m.orders(s.data, t.data, p.d.mat, gbsm.model="nb",
#'   metric="raw", orders = c(3:5, 8, 10, 15, 20), d.f=4,
#'    degree=3, n=1000, k=5, p=0.8, type="k-fold", scat.plots=TRUE,
#'     response.curves=TRUE, j.occs.distrbn=TRUE, mp.plots=TRUE,
#'      start=seq(-0.1, 0, length.out=(ncol(t.data)+2)*4+1))
#'
#'  jp4$contbn_table[[1]]
#'  jp4$model.validation.table
#'  jp4$jo.orders
#'
#' ## Close the open plots.gbsm.pdf file before running the 3rd example
#'  RNGkind(sample.kind = "Rejection")
#'  set.seed(1)
#'  jp5 <- msco::gbsm_m.orders(s.data, t.data, p.d.mat, gbsm.model="quasipoisson",
#'   metric="raw", orders = c(3:5, 8, 10, 15, 20), d.f=4,
#'    degree=3, n=1000, k=5, p=0.8, type="k-fold", scat.plots=TRUE,
#'     response.curves=TRUE, j.occs.distrbn=TRUE, mp.plots=TRUE,
#'      start=seq(-0.1, 0, length.out=(ncol(t.data)+2)*4+1))
#'
#'  jp5$contbn_table[[1]]
#'  jp5$model.validation.table
#'  jp5$jo.orders
#'
#'  }
#'
#' @export
#' @md

gbsm_m.orders <- function(s.data, t.data, p.d.mat, metric="Simpson_eqn", orders, d.f=4, degree=3, n=1000, k=5, p=0.8, type="k-fold", gbsm.model,
                          scat.plots=FALSE, response.curves=TRUE, j.occs.distrbn=FALSE, mp.plots=FALSE, start=seq(-0.1, 0, length.out=(ncol(t.data)+2)*4+1)){

  grDevices::pdf(file = paste0(system.file("ms", package = "msco"), "/plots.gbsm.pdf"), height = 8.27, width = 4)
  graphics::par(mar=c(5,5,4,1)+.1)
  graphics::par(mfrow=c((length(orders)+1)/2,2))

  vars2 <- rep(NA, length(orders))
  pred.variables <- c(names(t.data), "P.dist", "E.rate")
  cols <- c("red","blue","black","green", "orange", "brown", "purple", "yellow", "pink")
  names.var <- c()
  mss <- list()
  contbn_tablee <- list()
  cvalid <- list()
  cvalid_TEs <- matrix(NA, nrow = length(orders), ncol = 3)
  order.names <- c()
  N <- ncol(s.data)
  for (i in orders) {
    contbn_table <- `names<-`(data.frame(rep(NA, ncol(t.data)+2), rep(NA, ncol(t.data)+2), rep(NA,ncol(t.data)+2),
                                         rep(NA,ncol(t.data)+2)), c("predictor", "var.expld_M1", "var.expld_M2",
                                                                    "contribution"))
    pred.cont <- list()
    if((which(orders==i)%%2) == 0){
      mss[[i]] <- gbsm(s.data, t.data, p.d.mat, d.f=d.f, metric=metric, order.jo=i, degree=degree, n=n, b.plots=FALSE,
                       gbsm.model = gbsm.model, response.curves=response.curves, ylabel=FALSE, scat.plot=FALSE, leg = 0, start = start)
    }else{
      mss[[i]] <- gbsm(s.data, t.data, p.d.mat, d.f=d.f, metric=metric, order.jo=i, degree=degree, n=n, b.plots=FALSE,
                       gbsm.model = gbsm.model, response.curves=response.curves, ylabel=TRUE, scat.plot=FALSE, leg = 0, start = start)
    }
    bs_pred <- mss[[i]]$bs_pred
    j.occs <- mss[[i]]$j.occs
    gof <- mss[[i]]$var.expld

    ## Cross-validation
    if(type=="validation.set"){
      cvalid[[i]] <- cross_valid(mss[[i]], type, k=k, p=p)
    }else{
      cvalid[[i]] <- cross_valid(mss[[i]], type, k=k, p=p)$results[2:4]
    }


    ### Predictor contribution
    for(j in 1:(ncol(t.data)+2)){
      data <- bs_pred[, -which(names(`names<-`(bs_pred,gsub("[[:digit:]]", "", names(bs_pred)))) %in% c(unique(
        names(`names<-`(bs_pred,gsub("[[:digit:]]", "", names(bs_pred)))))[j]))]

      if((metric %in% c("raw_prop", "Simpson_eqn", "Sorensen_eqn"))==TRUE){
        pred.cont[[j]] <- suppressWarnings(glm2::glm2(j.occs ~ ., family=stats::quasibinomial(link="log"), data = data,
                                                      start = seq(-0.1, 0, length.out=ncol(data)+1)))
      }else if(metric=="raw" & (metric %in% c("raw_prop", "Simpson_eqn", "Sorensen_eqn"))!=TRUE & gbsm.model=="quasipoisson"){
        pred.cont[[j]] <- suppressWarnings(glm2::glm2(j.occs ~ ., family=stats::quasipoisson(link="log"), data = data,
                                                      start = seq(-0.1, 0, length.out=ncol(data)+1)))
      }else if(metric=="raw" & (metric %in% c("raw_prop", "Simpson_eqn", "Sorensen_eqn"))!=TRUE & gbsm.model=="nb"){
        pred.cont[[j]] <- suppressWarnings(MASS::glm.nb(j.occs ~ ., link = log, data = data, start = seq(-0.1, 0, length.out=ncol(data)+1)))
      }else  if(metric=="raw" & (metric %in% c("raw_prop", "Simpson_eqn", "Sorensen_eqn"))!=TRUE & (gbsm.model %in% c("quasipoisson", "nb"))!=TRUE){
        stop("Wrong gbsm.model used for 'raw' version of joint occupancy. It must either be 'quasipoisson' or 'nb'.")
      }

      contbn_table$predictor[j] <- unique(names(`names<-`(bs_pred,gsub("[[:digit:]]", "", names(bs_pred)))))[j]
      contbn_table$var.expld_M2[j] <- stats::cor(j.occs, as.numeric(suppressWarnings(stats::predict.glm(pred.cont[[j]],
                                                                                             newdata = data,
                                                                                             type = "response"))))^2
      contbn_table$var.expld_M1[j] <- gof
      contbn_table$contribution[j] <- (contbn_table$var.expld_M1[j] - contbn_table$var.expld_M2[j])/(contbn_table$var.expld_M1[j])
    }
    contbn_tablee[[i]] <- contbn_table
    vars2[i] <- gof

    order.names[i] <- paste("order", i)
  }

  order.names <- order.names[stats::complete.cases(order.names)]
  contbn_tablee <- contbn_tablee[!sapply(contbn_tablee,is.null)]
  contbn_tablee <- `names<-`(contbn_tablee, order.names)

  if(response.curves==TRUE){
    graphics::plot.new()
    nvar <- 1
    ncex <- 1
    if(length(pred.variables)>9){
      nvar <- 2
      ncex <- 0.6
    }
    graphics::legend("center", legend = pred.variables, col = cols, bty = "n", lty=1:length(pred.variables),
                     lwd=2, pch = 1:(length(t.data)+2), cex = ncex, ncol = nvar)
  }

  if(scat.plots==TRUE){
    graphics::par(mar=c(5,5,4,1)+.1)
    graphics::par(mfrow=c((length(orders)+1)/2,2))
    for (kv in orders) {
      if((which(orders==kv)%%2) == 0){
        plot(mss[[kv]]$j.occs, (mss[[kv]]$pred.j.occs), xlab="J. occupancy", ylab=" ",
             main = noquote(paste("Order", kv)))
      }else{
        plot(mss[[kv]]$j.occs, (mss[[kv]]$pred.j.occs), xlab="J. occupancy", ylab="Predicted J.occ",
             main = noquote(paste("Order", kv)))
      }
    }
  }
  if(j.occs.distrbn==TRUE){
    graphics::par(mar=c(5,5,4,1)+.1)
    graphics::par(mfrow=c((length(orders)+1)/2,2))
    occs <- rowSums(s.data)
    occs_scaled <- (occs-min(occs))/(max(occs)-min(occs)) ## scaled to be in [0,1]
    for (kl in c(1,orders)) {
      if(kl==1){
        graphics::hist(occs_scaled, xlab="Occupancy", ylab="Frequency", main = noquote(paste("Order", 1)),
                       xlim=c(0,1), breaks=seq(0, 1, 0.1))
      }else{
        if(metric!="raw"){
          if((which(orders==kl)%%2) == 0){
            graphics::hist(mss[[kl]]$j.occs, xlab="J. occupancy", ylab="Frequency", main = noquote(paste("Order", kl)),
                           xlim=c(0,1), breaks=seq(0, 1, 0.1))
          }else{
            graphics::hist(mss[[kl]]$j.occs, xlab="J. occupancy", ylab=" ", main = noquote(paste("Order", kl)),
                           xlim=c(0,1), breaks=seq(0, 1, 0.1))
          }
        }else if(metric=="raw"){
          if((which(orders==kl)%%2) == 0){
            graphics::hist((mss[[kl]]$j.occs-min(mss[[kl]]$j.occs))/(max(mss[[kl]]$j.occs)-min(mss[[kl]]$j.occs)),
                           xlab="J. occupancy", ylab="Frequency", main = noquote(paste("Order", kl)),
                           xlim=c(0,1), breaks=seq(0, 1, 0.1))
          }else{
            graphics::hist((mss[[kl]]$j.occs-min(mss[[kl]]$j.occs))/(max(mss[[kl]]$j.occs)-min(mss[[kl]]$j.occs)),
                           xlab="J. occupancy", ylab=" ", main = noquote(paste("Order", kl)),
                           xlim=c(0,1), breaks=seq(0, 1, 0.1))
          }

        }
      }
    }
  }

  if(length(which(is.na(vars2)==TRUE)==TRUE)>=1){
    vars2 <- vars2[stats::complete.cases(vars2)]
  }

  GoFs <- `names<-`(data.frame(vars2), "Rsquared_gf")


  for (j in 1:length(orders)) {
    cvalid_TEs[j,] <- as.numeric(Filter(Negate(is.null), cvalid)[[j]])
  }
  cvalid_rsq <- cvalid_TEs[,2]
  cvalid_rsq <- `names<-`(as.data.frame(cvalid_rsq), "Rsquared_cv")
  mod.valid <- cbind(orders, GoFs, cvalid_rsq)

  ### Model performance plots
  if(mp.plots==TRUE){
    GoFs <- GoFs$Rsquared_gf
    cvalid_rsq <- cvalid_rsq$Rsquared_cv
    ## Fit GoFs (with Orders) to exponential-power law
    mod2 <- minpack.lm::nlsLM(GoFs~a*exp(b*orders)*orders^c, start=list(a=1, b=0, c=0))
    pred.GoFs <- stats::predict(mod2, data.frame(orders), type="response")
    cc <- summary(mod2)$coefficients[,1]

    ## Fit model performance from cross validation (cvalid_rsq; with Orders) to exponential-power law
    mod3 <- minpack.lm::nlsLM(cvalid_rsq~a*exp(b*orders)*orders^c, start=list(a=1, b=0, c=0))
    pred.cvalid_rsq <- stats::predict(mod3, data.frame(orders), type="response")
    cc2 <- summary(mod3)$coefficients[,1]

    ## Compare GoFs and cvalid_rsq to have overall model performance
    mod4 <- stats::lm(GoFs~cvalid_rsq)
    mod4.val <- seq(from=range(cvalid_rsq)[1], to=range(cvalid_rsq)[2], length.out=1000)
    mod4y <- (summary(mod4)$coefficients[,1][[2]]*mod4.val) + summary(mod4)$coefficients[,1][[1]]

    ############ Transform using seq()
    orders.trans <- seq(from=range(orders)[1], to=range(orders)[2], length.out=1000)
    GoFs.est <- (cc[[1]])*(exp(cc[[2]]*orders.trans)) * (orders.trans)^(cc[[3]])
    cvalid_rsq.est <- (cc2[[1]])*(exp(cc2[[2]]*orders.trans)) * (orders.trans)^(cc2[[3]])

    cols <- c("red","blue","black","green", "orange", "brown", "purple")
    graphics::par(mar=c(5,5,4,1)+.1)
    graphics::par(mfrow=c(2,1))

    graphics::plot(orders.trans, GoFs.est, type = "l", lwd=2, lty=1, col=cols[1], main = "Effect of predictors on J. occ",
                   ylim=range(GoFs, cvalid_rsq, GoFs.est, cvalid_rsq.est),
                   xlab = "Joint occupancy order", ylab = noquote(expression(paste("Variance explained", "  ", (r^2)))))
    graphics::points(orders, GoFs, lwd=2, col=cols[1], pch=1)
    graphics::lines(orders.trans, cvalid_rsq.est, type="l", lty=2, lwd=2, col=cols[2])
    graphics::points(orders, cvalid_rsq, lwd=2, col=cols[2], pch=4)
    graphics::legend("topright", legend = c("Rsquared_gf", "Rsquared_cv"), pch=c(1,4), col = cols, lty=c(1,2), lwd=2,
                     bty = "n", cex=0.7, ncol = 1)
    graphics::mtext("(a)", side = 3, adj = -0.2, line = 1.5, cex=1.2, font = 2)

    graphics::plot(mod4.val, mod4y, type="l", col="black", lwd=2, main = "Model performance", ylim=range(cvalid_rsq, mod4y),
                   xlab = "Rsquared_cv", ylab = "Rsquared_gf")
    graphics::points(GoFs, cvalid_rsq, lwd=2, col=cols[3])
    graphics::mtext("(b)", side = 3, adj = -0.2, line = 1.5, cex = 1.2, font = 2)
    graphics::text(0.29, (max(mod4.val)-0.05), paste("mp", "", "=", "", round(((stats::cor(GoFs, cvalid_rsq))^2)*100, 1),"%"), font=2)

  }
  grDevices::dev.off()

  m.orders <- list()
  m.orders$contbn_table <- contbn_tablee
  m.orders$model.validation.table <- mod.valid
  m.orders$metric <- metric
  m.orders$d.f <- d.f
  m.orders$n <- n
  m.orders$degree <- degree
  m.orders$jo.orders <- orders
  m.orders$gbsm.plots <- print(noquote("Check msco's 'inst/ms' directory in your R library for a 'plots.gbsm.pdf' file."))
  return(m.orders)
}

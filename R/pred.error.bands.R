#' Prediction uncertainty
#'
#' This function plots the response curves showing the effect of the predictors (i.e. trait-based and
#'  neutral forces) on joint occupancy as the response variable, with prediction error bands (as the
#'   standard deviation from the mean of the response variable) for all orders of joint occupancy.
#'
#' @param s.data A species-by-site presence/absence `data.frame` with entries indicating
#' occurrence (1) and non-occurrence (0) of species in a site.
#' @param t.data A `data.frame` with traits as columns and species as rows. The species must be the same as
#'  in `s.data`.
#' @param p.d.mat A symmetric `matrix` with `dimnames` as species and entries indicating the
#'  phylogenetic distance between any two of them (species).
#' @param orders As for \link[msco]{gbsm_m.orders}
#' @param degree As for \link[msco]{gbsm_m.orders}.
#' @param start As for \link[msco]{gbsm_m.orders}.
#' @param n As for \link[msco]{gbsm_m.orders}.
#' @param d.f As for \link[msco]{gbsm_m.orders}.
#' @param metric As for \link[msco]{gbsm_m.orders}.
#' @param simm Number of Monte Carlo simulations performed
#'
#' @return `pred.error.bands` function returns:
#' \item{`predictors`}{a `data.frame` of predictors}
#' \item{`responses`}{a `data.frame` of response values of predictors}
#' \item{`responses.sim_stats`}{a `data.frame` of the reponses' mean and standard deviation
#'  (from `simm` replicates), and}
#' * the response curves with prediction error bands for all orders of joint occupancy
#'
#' @references
#' \enumerate{
#'  \item{Lagat, V. K., Latombe, G. and Hui, C. (2021a). *A multi-species co-occurrence
#'  index to avoid type II errors in null model testing*. DOI: `<To be added>`.}
#'
#'  \item{Lagat, V. K., Latombe, G. and Hui, C. (2021b). *Dissecting the effects of
#'   neutral encounter versus functional traits on multi-order species interactions
#'    and co-occurrence with generalised B-spline modelling*. DOI: `<To be added>`.}
#' }
#' @examples
#' \dontrun{
#'  my.path <- system.file("extdata/gsmdat", package = "msco")
#'  setwd(my.path)
#'  s.data <- get(load("s.data.csv")) ## Species-by-site matrix
#'  t.data <- get(load("t.data.csv")) ## Species-by-Trait matrix
#'  p.d.mat <- get(load("p.d.mat.csv")) ## Species-by-species phylogenetic distance matrix
#'
#'  RNGkind(sample.kind = "Rejection")
#'  set.seed(0)
#'  pe <- msco::pred.error.bands(s.data, t.data, p.d.mat, metric="Simpson_eqn", d.f=4, simm=10,
#'   orders = c(2:5, 8, 10, 15), degree=3, n=1000,
#'    start=seq(-0.1, 0, length.out=(ncol(t.data)+2)*4+1))
#'
#'  pe$predictors$`order 2`
#'  pe$responses$`order 2`
#'  pe$responses.sim_stats$`order 2`
#'
#'  pe$predictors$`order 3`
#'  pe$responses$`order 3`
#'  pe$responses.sim_stats$`order 3`
#'
#'  pe$predictors$`order 10`
#'  pe$responses$`order 10`
#'  pe$responses.sim_stats$`order 10`
#'
#'  }
#'
#' @export
#' @md

pred.error.bands <- function(s.data, t.data, p.d.mat, metric="Simpson_eqn", d.f=4, simm=10, orders, degree=3, n=1000, start=seq(-0.1, 0, length.out=(ncol(t.data)+2)*4+1)){

  grDevices::pdf(file = paste0(system.file("ms", package = "msco"), "/pred.error.bands.pdf"), paper="a4r", height = 8.27, width = 11.69)
  graphics::par(mar=c(4,4,2,0.5)+.1)
  graphics::par(mfcol=c((ncol(t.data)+2),length(orders)))

  Predictors <- list()
  Responses <- list()
  Responses_stats <- list()
  order.names <- c()
  for (i in orders) {
    if(i==2){
      pr <- pred.error.bands.table(s.data, t.data, p.d.mat, metric, d.f, simm, order.jo=i, degree, n, start)
      predictors <- pr$predictors
      responses <- pr$responses
      response_stats <- pr$respns_dispn.table

      order.names[i] <- paste("order", i)

      Predictors[[i]] <- predictors
      Responses[[i]] <- responses
      Responses_stats[[i]] <- response_stats


      graphics::plot(predictors[,1], response_stats[, (1+(1-1)*2)], type= "l", main=paste("Order", i), lwd=2, ylab = expression(bold("Effect on J. Occ")), xlab = noquote(names(predictors)[1]),
                     ylim=range((response_stats[,(1+(1-1)*2)] - response_stats[,(2+(1-1)*2)]), (response_stats[,(1+(1-1)*2)] + response_stats[,(2+(1-1)*2)])))

      graphics::polygon(c(rev(predictors[,1]), predictors[,1]), c(rev((response_stats[,(1+(1-1)*2)] - response_stats[,(2+(1-1)*2)])),
                                                                  (response_stats[,(1+(1-1)*2)] + response_stats[,(2+(1-1)*2)])), col = 'grey80', border = NA)
      graphics::lines(predictors[,1], response_stats[, (1+(1-1)*2)], type= "l", lwd=2)

      for (v in 2:ncol(predictors)) {
        graphics::plot(predictors[,v], response_stats[, (1+(v-1)*2)], type= "l", lwd=2, ylab = expression(bold("Effect on J. Occ")), xlab = noquote(names(predictors)[v]),
                       ylim=range((response_stats[,(1+(v-1)*2)] - response_stats[,(2+(v-1)*2)]), (response_stats[,(1+(v-1)*2)] + response_stats[,(2+(v-1)*2)])))

        graphics::polygon(c(rev(predictors[,v]), predictors[,v]), c(rev((response_stats[,(1+(v-1)*2)] - response_stats[,(2+(v-1)*2)])),
                                                                    (response_stats[,(1+(v-1)*2)] + response_stats[,(2+(v-1)*2)])), col = 'grey80', border = NA)
        graphics::lines(predictors[,v], response_stats[, (1+(v-1)*2)], type= "l", lwd=2)
      }
    }else{
      pr <- pred.error.bands.table(s.data, t.data, p.d.mat, metric, d.f, simm, order.jo=i, degree, n, start)
      responses <- pr$responses
      predictors <- pr$predictors
      response_stats <- pr$respns_dispn.table

      order.names[i] <- paste("order", i)

      Predictors[[i]] <- predictors
      Responses[[i]] <- responses
      Responses_stats[[i]] <- response_stats

      graphics::plot(predictors[,1], response_stats[, (1+(1-1)*2)], type= "l", main=paste("Order", i), lwd=2, ylab = "", xlab = noquote(names(predictors)[1]),
                     ylim=range((response_stats[,(1+(1-1)*2)] - response_stats[,(2+(1-1)*2)]), (response_stats[,(1+(1-1)*2)] + response_stats[,(2+(1-1)*2)])))

      graphics::polygon(c(rev(predictors[,1]), predictors[,1]), c(rev((response_stats[,(1+(1-1)*2)] - response_stats[,(2+(1-1)*2)])),
                                                                  (response_stats[,(1+(1-1)*2)] + response_stats[,(2+(1-1)*2)])), col = 'grey80', border = NA)
      graphics::lines(predictors[,1], response_stats[, (1+(1-1)*2)], type= "l", lwd=2)

      for (v in 2:ncol(predictors)) {
        graphics::plot(predictors[,v], response_stats[, (1+(v-1)*2)], type= "l", lwd=2, ylab = "", xlab = noquote(names(predictors)[v]),
                       ylim=range((response_stats[,(1+(v-1)*2)] - response_stats[,(2+(v-1)*2)]), (response_stats[,(1+(v-1)*2)] + response_stats[,(2+(v-1)*2)])))

        graphics::polygon(c(rev(predictors[,v]), predictors[,v]), c(rev((response_stats[,(1+(v-1)*2)] - response_stats[,(2+(v-1)*2)])),
                                                                    (response_stats[,(1+(v-1)*2)] + response_stats[,(2+(v-1)*2)])), col = 'grey80', border = NA)
        graphics::lines(predictors[,v], response_stats[, (1+(v-1)*2)], type= "l", lwd=2)
      }
    }
  }
  grDevices::dev.off()
  base::system(paste0('open "', paste0(system.file("ms", package = "msco"), "/pred.error.bands.pdf"), '"'))

  order.names <- order.names[stats::complete.cases(order.names)]
  Predictors <- Predictors[!sapply(Predictors,is.null)]
  Responses <- Responses[!sapply(Responses,is.null)]
  Responses_stats <- Responses_stats[!sapply(Responses_stats,is.null)]

  peb <- list()
  peb$predictors <- `names<-`(Predictors, order.names)
  peb$responses <- `names<-`(Responses, order.names)
  peb$responses.sim_stats <- `names<-`(Responses_stats, order.names)
  return(peb)
}



pred.error.bands.table <- function(s.data, t.data, p.d.mat, metric="Simpson_eqn", d.f=4, simm=50, order.jo=2,
                                   degree=3, n=1000, max.vif=40, max.vif2=30, start.range=c(-0.1,0)){

  predictors <- msco::gbsm(s.data, t.data, p.d.mat, metric = metric, d.f=d.f, order.jo=order.jo, max.vif = max.vif, max.vif2 = max.vif2,
                           degree=degree, n=n, b.plots=FALSE, scat.plot=FALSE, response.curves=FALSE, leg=1, start.range=start.range)$Predictors

  simm_respns.table <- `names<-`(as.data.frame(matrix(NA, nrow = nrow(predictors), ncol = (ncol(t.data)+2)*simm)), rep(c(names(t.data), "P.dist", "E.rate"),simm))

  for (ki in 1:simm) {
    responses <- msco::gbsm(s.data, t.data, p.d.mat, metric = metric, d.f=d.f, order.jo=order.jo, max.vif = max.vif, max.vif2 = max.vif2,
                            degree=degree, n=n, b.plots=FALSE, scat.plot=FALSE, response.curves=FALSE, leg=1, start.range=start.range)$Responses

    simm_respns.table[,(1+(ki-1)*5):(5+(ki-1)*5)] <- responses

  }

  responses_mean <- c()
  responses_sd <- c()
  respns_dispn.table <- as.data.frame(matrix(NA, nrow = nrow(predictors), ncol=ncol(predictors)*2))
  for (i in 1:ncol(predictors)) {
    responses_sd <- apply(simm_respns.table[,which(names(simm_respns.table)==names(simm_respns.table)[i])],1, stats::sd)
    responses_mean <- apply(simm_respns.table[,which(names(simm_respns.table)==names(simm_respns.table)[i])],1, mean)
    respns_dispn.table[,(1+(i-1)*2):(2+(i-1)*2)] <- cbind(responses_mean, responses_sd)

  }
  mynames <- c()
  for (i in 1:ncol(predictors)) {
    mynames[(1+(i-1)*2):(2+(i-1)*2)] <-  c(paste0(names(predictors)[i], "_", "mean"), paste0(names(predictors)[i], "_", "sd"))
  }
  respns_dispn.table <- `names<-`(respns_dispn.table, mynames)

  prr <- list()
  prr$respns_dispn.table <- respns_dispn.table
  prr$predictors <- predictors
  prr$responses <- responses
  return(prr)
}


mat.reshfler <- function(Jo.eng_Obj){

  simdata <- Jo.eng_Obj$Sim.Data
  s.data <- Jo.eng_Obj$Obs.Data
  if(is.null(s.data)){
    stop("Your 'Obs.data' in 'Jo.eng_Obj' is set to FALSE. Set it to TRUE then try again.")
  }
  if(is.null(simdata)){
    stop("Your 'Sim.Data' in 'Jo.eng_Obj' is set to FALSE. Set it to TRUE then try again.")
  }

  graphics::par(mfrow = c(1,2))

  ## OBSERVED

  plot(s.data, xlim = c(0,ncol(s.data)), ylim = c(0, nrow(s.data)), type = "n", ann = FALSE, axes = FALSE)
  graphics::mtext("Sites", side = 1, font = 2)
  graphics::mtext("Species", side = 2, font = 2)
  graphics::mtext("Observed", side = 3, font = 2, col = "grey30")
  yr <- rep(0:(nrow(s.data)-1), ncol(s.data))
  xr <- rep(0:(ncol(s.data)-1), each = nrow(s.data))
  box.cols <- c("white","grey30")
  Col_Vec <- box.cols[as.integer(s.data)+1]
  graphics::rect(xr, yr, xr+1, yr+1, col = Col_Vec)

  ## SIMULATED

  plot(simdata, xlim = c(0, ncol(simdata)), ylim = c(0,nrow(simdata)), type = "n", ann = FALSE, axes = FALSE)
  graphics::mtext("Sites", side = 1, font=2)
  graphics::mtext("Species", side = 2, font = 2)
  graphics::mtext("Simulated", side = 3, font = 2, col = "darkgreen")
  yr <- rep(0:(nrow(simdata)-1), ncol(simdata))
  xr <- rep(0:(ncol(simdata)-1), each = nrow(simdata))
  box.cols <- c("white","darkgreen")
  Col_Vec <- box.cols[as.integer(simdata)+1]
  graphics::rect(xr, yr, xr+1, yr+1, col = Col_Vec)
}




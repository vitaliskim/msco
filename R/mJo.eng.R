#' Joint occupancy model engine for multiple communities
#'
#' This function is the engine behind the null model testing of species co-occurrence patterns, and
#'  analyses of the joint occupancy decline and the parametric forms of this decline, for multiple
#'   communities. In particular:
#' \itemize{
#' \item{It performs the null model testing of species co-occurrence patterns and generates the
#'  archetypes depicting how joint occupancy declines with the number of species (the order
#'   of msco) based on species-by-site presence/absence `.csv` data matrices. From these archetypes,
#'    inferences can be made according to the implemented null models};
#'  \item{Determines the robustness of the exponential, power law and exponential-power law forms of
#'   joint occupancy decline by computing the Pearson's \eqn{r^2} between the joint occupancy values
#'    of the observed data and predicted data, for all orders of species};
#'  \item{Gives a summary of the total number of communities (under each and for all archetypes) whose
#'   forms of joint occupancy decline have \eqn{r^2 > 0.95}};
#'  \item{Computes the AIC and Delta AIC of joint occupancy decline regression models for all communities};
#'  \item{Computes the total number of communities:
#'    \itemize{
#'   \item{with exponential as the best form of joint occupancy decline than power law and vice versa;}
#'   \item{with either of the three regression models (exponential, power law and exponential-power law)
#'    having the best form of the joint occupancy decline;}}}
#'  \item{Estimates the parameters of:
#'   \enumerate{
#'    \item{**exponential:** **\eqn{j^{\{i\}} = a \times exp(b \times i)}**;}
#'    \item{**power law:** **\eqn{j^{\{i\}} = a \times i^b}**; and}
#'    \item{**exponential-power law:** **\eqn{j^{\{i\}} = a \times exp(b \times i) \times i^c}**}
#'  }}forms of joint occupancy decline, respectively, and their 95% confidence interval.}
#'
#'@details `mJo.eng` function is useful when analyzing multiple species-by-site presence/absence
#'  data matrices at once. If one community matrix is analyzed, the outputs of the function
#'   \link[msco]{Jo.eng} should suffice.
#'
#' @param my.files A vector containing names of species-by-site presence/absence `.csv` data matrices.
#'  The data matrices should be saved in the working directory.
#' @param algo Simulation algorithm used. The possible options to choose from are: `sim1`,
#'  `sim2`, `sim3`, `sim4`, `sim5`, `sim6`, `sim7`, `sim8`, and `sim9`, all from
#'   Gotelli (2000). `sim2` is highly recommended (see Lagat *et al.,* 2021a).
#' @param metric The type of rescaling applied to the joint occupancy metric. Available options are:
#'  `Simpson_eqn` for Simpson equivalent, `Sorensen_eqn` for Sorensen equivalent, and `raw` for the
#'   raw form of index without rescaling.
#' @param nReps Number of simulations used in the null model test.
#' @param Archetypes A Boolean indicating if the archetypes of the patterns
#'  of species co-occurrences in multiple communities should be included in the output.
#' @param AICs A Boolean indicating whether the akaike information criterion (AIC) and Delta AIC
#'  of joint occupancy decline regression models for all communities should be included in the output.
#' @param params A Boolean indicating whether parameter estimates of the joint occupancy decline
#'  regression models should be included in the output.
#' @param best.mod2 A Boolean indicating if exponential and power law regression model comparisons
#'  should be included in the output.
#' @param best.mod3 A Boolean indicating if exponential, power law and exponential-power law
#'  regression model comparisons should be included in the output.
#' @param params_c.i A Boolean indicating if 95% C.I of the parameter estimates of the joint occupancy
#'  decline regression models should be included in the output.
#' @param my.r2 A Boolean indicating if the robustness of joint occupancy decline regression models
#'  should be computed and output.
#' @param my.r2.s A Boolean indicating if the robustness summary values of joint occupancy decline
#'  regression models should be computed and output.
#' @return `mJo.eng` function returns a list containing the following outputs:
#'
#' $`Archs`
#'
#' For every community, a `list` consisting of:
#'
#' \itemize{
#'   \item `$nmod_stats`: A data frame with the summary statistics for the null model test; and
#'   \item `$Archetype`: Archetypes of the patterns of species co-occurrences in ecological
#'    communities/matrices (`my.files`). These archetypes must be \eqn{\in \{}"A1",
#'     "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9"\eqn{\}} or "NA". "NA" could be the
#'      combinations of two or more of the nine expected archetypes.}
#' $`all.AICs`
#'
#' A `list` of `data.frame`s containig the following components:
#' \item{`df`}{The number of parameters in each of the three (exponential, power law and
#'     exponential-power law) joint occupancy decline regression models.}
#' \item{`aic`}{The aic values for each of the three joint occupancy decline regression models.}
#' \item{`delta_aic3`}{The `delta_aic` values for each of the three joint occupancy decline regression
#'      models.}
#' \item{`delta_aic2`}{The `delta_aic` values for exponential and power law forms of joint occupancy
#'     decline regression models.}
#' $`params`
#'
#' A `data.frame` consisting of:
#' \item{`arch`}{The archetypes of the patterns of species co-occurrences in each of the  species-by-site
#'     presence/absence `.csv` data matrices.}
#' \item{`a.ex`}{The `a` parameter estimate of the exponential form of joint occupancy decline.}
#' \item{`b.ex`}{The `b` parameter estimate of the exponential form of joint occupancy decline.}
#' \item{`a.pl`}{The `a` parameter estimate of the power law form of joint occupancy decline.}
#' \item{`b.pl`}{The `b` parameter estimate of the power law form of joint occupancy decline.}
#' \item{`a.expl`}{The `a` parameter estimate of the exponential-power law form of joint occupancy decline.}
#' \item{`b.expl`}{The `b` parameter estimate of the exponential-power law form of joint occupancy decline.}
#' \item{`c.expl`}{The `c` parameter estimate of the exponential-power law form of joint occupancy decline.}
#' $`best.mod2`
#'
#' A`table` containig the following components:
#' \item{`n`}{The number of ecological communities represented by species-by-site
#'      presence/absence `.csv` data matrices.}
#' \item{`n.lwst_aic`}{The number of communities with exponential as the best
#'     form of joint occupancy decline than power law.}
#' \item{`n.delta_aic`}{The number of communities whose exponential and power law forms of joint occupancy
#'      decline have `delta_aic = 0`, respectively. This number must be equal to `n.lwst_aic`.}
#' \item{`%`}{The percentage of `n.lwst_aic` (or `n.delta_aic`) relative to the total number of
#'      communities (`n`) analyzed.}
#' $`best.mod3`
#'
#' A `table` containig the following components:
#' \item{`n`}{The number of ecological communities represented by species-by-site
#'     presence/absence `.csv` data matrices.}
#' \item{`n.lwst_aic`}{The number of communities with exponential or power law or exponential-power
#'  law as the best form of joint occupancy decline among the three (exponential, power law and
#'   exponential-power law) regression models.}
#' \item{`n.delta_aic`}{The number of communities whose exponential, power law and exponential-power
#'     law forms of joint occupancy decline, respectively, have `delta_aic = 0`. This number must be
#'    equal to `n.lwst_aic`.}
#' \item{`%`}{The percentage of `n.lwst_aic` (or `n.delta_aic`) relative to the total number of
#'    communities (`n`) analyzed.}
#' $`params_c.i`
#'
#' A `data.frame` consisting of:
#' \item{`arch`}{The archetypes of the patterns of species co-occurrences in each of the  species-by-site
#'     presence/absence `.csv` data matrices.}
#' \item{`n`}{The number of communities under every archetype.}
#' \item{`ex_%`}{The percentages of the number of communities (under every archetype) where
#'     exponential form of joint occupancy decline fitted better than power law.}
#' \item{`a.ex`}{The 95% closed confidence interval of the `a` parameter estimates of the exponential
#'     form of joint occupancy decline, under every archetype.}
#' \item{`b.ex`}{The 95% closed confidence interval of the `b` parameter estimates of the exponential
#'     form of joint occupancy decline, under every archetype.}
#' \item{`p.l_%`}{The percentages of the number of communities (under every archetype) where
#'     power law form of joint occupancy decline fitted better than exponential.}
#' \item{`a.pl`}{The 95% closed confidence interval of the `a` parameter estimates of the power law
#'     form of joint occupancy decline, under every archetype.}
#' \item{`b.pl`}{The 95% closed confidence interval of the `b` parameter estimates of the power law
#'      form of joint occupancy decline, under every archetype.}
#' \item{`ex.pl_%`}{The percentages of the number of communities (under every archetype) where exponential-power
#'      law form of joint occupancy decline fitted better than both the exponential and power law forms.}
#' \item{`a.expl`}{The 95% closed confidence interval of the `a` parameter estimates of the exponential-power law
#'     form of joint occupancy decline, under every archetype.}
#' \item{`b.expl`}{The 95% closed confidence interval of the `b` parameter estimates of the exponential-power law
#'     form of joint occupancy decline, under every archetype.}
#' \item{`c.expl`}{The 95% closed confidence interval of the `c` parameter estimates of the exponential-power law
#'     form of joint occupancy decline, under every archetype.}
#' $`r2`
#'
#' A `list` of `data.frame`s containig the following components:
#' \item{`rsq.ex`}{\eqn{r^2} for the exponential form of joint occupancy decline.}
#' \item{`rsq.pl`}{\eqn{r^2} for the power law form of joint occupancy decline.}
#' \item{`rsq.ex.pl`}{\eqn{r^2} for the exponential-power law form of joint occupancy decline.}
#' $`r2.s`
#' * A `list` containig the following components:
#'
#'   $`rsq.per.Archs`
#'
#'   + `Archs`: Archetypes of the patterns of species co-occurrences in each
#'     of the species-by-site presence/absence .csv data matrices.
#'   + `n.a`: Number of communities under each archetype.
#'   + `rsq.ex`: Number of communities under each archetype whose exponential
#'     forms of joint occupancy decline have \eqn{r^2 > 0.95}.
#'   + `rsq.pl`: Number of communities under each archetype whose power
#'    law forms of joint occupancy decline have \eqn{r^2 > 0.95}.
#'   + `rsq.ex-pl`: Number of communities under each archetype whose
#'     exponential-power law forms of joint occupancy decline have \eqn{r^2 > 0.95}.
#'
#'   $`rsq.all.Communities`
#'
#'   + `n`: Number of all communities analyzed
#'   + `ex`: Number of communities whose exponential forms of joint occupancy
#'     decline have \eqn{r^2 > 0.95}
#'   + `pl`: Number of communities whose power law forms of joint occupancy
#'    decline have \eqn{r^2 > 0.95}
#'   + `ex.pl`: Number of communities whose exponential-power law forms
#'     of joint occupancy decline have \eqn{r^2 > 0.95}
#'
#' $`m.Jo.plots`
#'
#' Produces a `.pdf` file with multiple figures each consisting of the following plots:
#'
#' \item{(a)}{as for \link[msco]{Jo.plots}}
#' \item{(b)}{as for \link[msco]{Jo.plots}}
#' \item{(c)}{as for \link[msco]{Jo.plots}}
#' \item{(d)}{as for \link[msco]{Jo.plots}}
#' \item{(e)}{as for \link[msco]{Jo.plots}}
#' @references
#' \enumerate{
#'  \item{Lagat, V. K., Latombe, G. and Hui, C. (2021a). *A multi-species co-occurrence
#'  index to avoid type II errors in null model testing*. DOI: `<To be added>`.}
#'
#'  \item{Gotelli, N. J. (2000). Null model analysis of species co-occurrence patterns.
#'  *Ecology, 81(9)*, 2606-2621. <https://doi.org/10.1890/0012-9658(2000)081[2606:NMAOSC]2.0.CO;2>}
#'
#'  \item{Pearson, K. (1895) VII. Note on regression and inheritance in the
#'  case of two parents. *proceedings of the royal society of London,* **58**:240-242.
#'   <https://doi.org/10.1098/rspl.1895.0041>}
#'
#'  \item{Petrossian, G.A., Maxfield, M (2018). An information theory approach to hypothesis testing in
#'   criminological research. *crime science,* 7(1), 2. <https://doi.org/10.1186/s40163-018-0077-5>}
#'  }
#' @examples
#' \dontrun{
#'
#' my.path <- system.file("extdata", package = "msco")
#' setwd(my.path)
#' my.files <- gtools::mixedsort(list.files(path = my.path, pattern = "*.csv"))
#' my.res <- msco::mJo.eng(my.files = my.files, algo = "sim2", Archetypes = TRUE,
#'              metric = "raw", nReps = 999, AICs = FALSE, params = FALSE,
#'              best.mod2 = FALSE, best.mod3 = FALSE, params_c.i = FALSE,
#'              my.r2 = FALSE, my.r2.s = FALSE)
#' my.res$Archs$`252.csv`
#'
#' my.path2 <- system.file("extdata/myCSVs", package = "msco")
#' setwd(my.path2)
#' my.files2 <- gtools::mixedsort(list.files(path = my.path2, pattern = "*.csv"))
#' my.res2 <- msco::mJo.eng(my.files = my.files2[250:255], algo = "sim2", Archetypes = FALSE,
#'               metric = "raw", nReps = 999, AICs = FALSE, params = TRUE,
#'               best.mod2 = FALSE, best.mod3 = FALSE, params_c.i = FALSE,
#'               my.r2 = FALSE, my.r2.s = FALSE)
#' my.res2
#'
#' my.path2 <- system.file("extdata/myCSVs", package = "msco")
#' setwd(my.path2)
#' my.files2 <- gtools::mixedsort(list.files(path = my.path2, pattern = "*.csv"))
#' my.res3 <- msco::mJo.eng(my.files = my.files2[250:255], algo = "sim2", Archetypes = FALSE,
#'               metric = "raw", nReps = 999, AICs = FALSE, params = FALSE,
#'               best.mod2 = FALSE, best.mod3 = FALSE, params_c.i = TRUE,
#'               my.r2 = FALSE, my.r2.s = FALSE)
#' my.res3
#'  }
#' @export
#' @md

mJo.eng <- function(my.files,
                      algo = "sim2",
                      metric = "raw",
                      nReps = 999,
                      Archetypes = FALSE,
                      AICs = FALSE,
                      params = FALSE,
                      best.mod2 = FALSE,
                      best.mod3 = FALSE,
                      params_c.i = FALSE,
                      my.r2 = FALSE,
                      my.r2.s = FALSE){

  if(length(gtools::mixedsort(list.files(path = getwd(), pattern = "*.csv")))==0){
    stop("No \`.csv\` binary matrices in your working directory. The \"my.files\" file path
         and the working directory should be the same.")
  }

  coe <- NULL
  Archs <- list()
  nmstats <- list()
  nm_arch <- list()
  all.AICs <- list()
  r2 <- list()
  myfiles = lapply(my.files, utils::read.csv, header=T)
  param <- matrix(NA, ncol = 8, nrow = length(myfiles))


  grDevices::pdf(file = paste0(system.file("ms", package = "msco"), "/mJo.plots.pdf"), paper="a4r", height = 8.27, width = 11.69)

  for (j in 1:length(myfiles)) {
    coe <- Jo.eng(myfiles[[j]], algo = algo, nReps = nReps, metric = metric)
    Archs[[j]] <- coe$Archetype ### Archetypes
    nmstats[[j]] <- coe$nmod_stats ### StatisticsTable
    nm_arch[[j]] <- `names<-`(list(nmstats[[j]], Archs[[j]]), c("nmod_stats", "Archetype"))
    all.AICs[[j]] <- coe$AIC ### all. AICs
    r2[[j]] <- coe$r2 ### rsq


    pagenum::setPagenum(my.files[1])
    print(coe$all.plots)
    pagenum::pagenum(num = my.files[j], text = "File: ", just = c("left", "bottom"))


    param[j,] <- c(noquote(coe$Archetype),
                   round(as.numeric(cbind(coe$jo.coeff[1,1:2],
                                          coe$jo.coeff[2,1:2],
                                          coe$jo.coeff[3,1:3])), 4)) ### paramss
  }
  grDevices::dev.off()

  for(vee in 1:length(r2)){
    r2[[vee]] <- `rownames<-`(`colnames<-`(as.data.frame(
      matrix(r2[[vee]], ncol = 3)), c("rsq.ex", "rsq.pl", "rsq.ex.pl")), " ")
  }

  paramss <- as.data.frame(matrix(param, nrow = nrow(param), ncol = ncol(param)))
  names(paramss) <- c("Arch","a.ex","b.ex","a.pl","b.pl","a.expl","b.expl","c.expl")

  ###############################################################################
  ################ Exponential and power law comparison (best.mod2) #############
  ###############################################################################

  all2.ex <- c()
  all2.pl <- c()
  all2.exd <- c()
  all2.pld <- c()
  for (m in 1:length(all.AICs)) {
    if(all.AICs[[m]][,2][1]<all.AICs[[m]][,2][2]){
      all2.ex[m] <- "ex"
    }else if(all.AICs[[m]][,2][2]<all.AICs[[m]][,2][1]){
      all2.pl[m] <- "pl"
    }

    ########### Delta_AIC #################

    if(all.AICs[[m]][,4][1]<all.AICs[[m]][,4][2]){
      all2.exd[m] <- "exd"
    }else if(all.AICs[[m]][,4][2]<all.AICs[[m]][,4][1]){
      all2.pld[m] <- "pld"
    }
  }

  n <- length(my.files)
  expo <- c(n,length(which(all2.exd=="exd")), length(which(all2.ex=="ex")), ((length(which(
    all2.ex=="ex")))/n)*100)
  names(expo) <- c("n", "n.Delta_AIC", "n.Lwst_AIC", "%")

  p.law <- c(n,length(which(all2.pld=="pld")), length(which(all2.pl=="pl")), ((length(which(
    all2.pl=="pl")))/n)*100)
  names(p.law) <- c("n", "n.Delta_AIC", "n.Lwst_AIC", "%")

  both.models <- rbind(expo, p.law)
  both.models <- round(as.table(both.models))

  ###########################################################################################
  ########## Exponential, power law and exponential-p.law comparisons (best.mod3) ###########
  ###########################################################################################

  all3.ex <- c()
  all3.pl <- c()
  all3.ex.pl <- c()
  all3.exd <- c()
  all3.pld <- c()
  all3.ex.pld <- c()
  for (p in 1:length(all.AICs)) {
    if(min(all.AICs[[p]][,2][1],all.AICs[[p]][,2][2],
           all.AICs[[p]][,2][3])==all.AICs[[p]][,2][1]){
      all3.ex[p] <- "ex3"
    }
  }
  for (q in 1:length(all.AICs)) {
    if(min(all.AICs[[q]][,2][1],all.AICs[[q]][,2][2],
           all.AICs[[q]][,2][3])==all.AICs[[q]][,2][2]){
      all3.pl[q] <- "pl3"
    }
  }
  for (r in 1:length(all.AICs)) {
    if(min(all.AICs[[r]][,2][1],all.AICs[[r]][,2][2],
           all.AICs[[r]][,2][3])==all.AICs[[r]][,2][3]){
      all3.ex.pl[r] <- "ex.pl3"
    }
  }

  ######### Delta_AIC ###################

  for (pd in 1:length(all.AICs)) {
    if(min(all.AICs[[pd]][,3][1],all.AICs[[pd]][,3][2],
           all.AICs[[pd]][,3][3])==all.AICs[[pd]][,3][1]){
      all3.exd[pd] <- "ex3d"
    }
  }
  for (qd in 1:length(all.AICs)) {
    if(min(all.AICs[[qd]][,3][1],all.AICs[[qd]][,3][2],
           all.AICs[[qd]][,3][3])==all.AICs[[qd]][,3][2]){
      all3.pld[qd] <- "pl3d"
    }
  }
  for (rd in 1:length(all.AICs)) {
    if(min(all.AICs[[rd]][,3][1],all.AICs[[rd]][,3][2],
           all.AICs[[rd]][,3][3])==all.AICs[[rd]][,3][3]){
      all3.ex.pld[rd] <- "ex.pl3d"
    }
  }

  expo3 <- c(n, length(which(all3.exd=="ex3d")), length(which(all3.ex=="ex3")), ((length(which(
    all3.ex=="ex3")))/n)*100)
  names(expo3) <- c("n", "n.Delta_AIC", "n.Lwst_AIC", "%")

  p.law3 <- c(n, length(which(all3.pld=="pl3d")), length(which(all3.pl=="pl3")), ((length(which(
    all3.pl=="pl3")))/n)*100)
  names(p.law3) <- c("n", "n.Delta_AIC", "n.Lwst_AIC", "%")

  exp_p.law <- c(n, length(which(all3.ex.pld=="ex.pl3d")), length(which(all3.ex.pl=="ex.pl3")),
                 ((length(which(all3.ex.pl=="ex.pl3")))/n)*100)
  names(exp_p.law) <- c("n", "n.Delta_AIC", "n.Lwst_AIC", "%")

  three.models <- rbind(expo3, p.law3, exp_p.law)
  three.models <- round(as.table(three.models))

  ###############################################################################################
  ############################ params_c.i #######################################################
  ###############################################################################################

  ### 95% closed confidence interval for all the parameter estimates ####

  A1a.ex <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,2][which(paramss[,1]=="A1")]), probs = c(0.025,0.975))),nrow = 1),2)
  A1b.ex <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,3][which(paramss[,1]=="A1")]), probs = c(0.025,0.975))),nrow = 1),2)
  A1a.pl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,4][which(paramss[,1]=="A1")]), probs = c(0.025,0.975))),nrow = 1),2)
  A1b.pl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,5][which(paramss[,1]=="A1")]), probs = c(0.025,0.975))),nrow = 1),2)
  A1a.expl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,6][which(paramss[,1]=="A1")]), probs = c(0.025,0.975))),nrow = 1),2)
  A1b.expl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,7][which(paramss[,1]=="A1")]), probs = c(0.025,0.975))),nrow = 1),2)
  A1c.expl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,8][which(paramss[,1]=="A1")]), probs = c(0.025,0.975))),nrow = 1),2)

  A2a.ex <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,2][which(paramss[,1]=="A2")]), probs = c(0.025,0.975))),nrow = 1),2)
  A2b.ex <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,3][which(paramss[,1]=="A2")]), probs = c(0.025,0.975))),nrow = 1),2)
  A2a.pl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,4][which(paramss[,1]=="A2")]), probs = c(0.025,0.975))),nrow = 1),2)
  A2b.pl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,5][which(paramss[,1]=="A2")]), probs = c(0.025,0.975))),nrow = 1),2)
  A2a.expl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,6][which(paramss[,1]=="A2")]), probs = c(0.025,0.975))),nrow = 1),2)
  A2b.expl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,7][which(paramss[,1]=="A2")]), probs = c(0.025,0.975))),nrow = 1),2)
  A2c.expl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,8][which(paramss[,1]=="A2")]), probs = c(0.025,0.975))),nrow = 1),2)

  A3a.ex <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,2][which(paramss[,1]=="A3")]), probs = c(0.025,0.975))),nrow = 1),2)
  A3b.ex <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,3][which(paramss[,1]=="A3")]), probs = c(0.025,0.975))),nrow = 1),2)
  A3a.pl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,4][which(paramss[,1]=="A3")]), probs = c(0.025,0.975))),nrow = 1),2)
  A3b.pl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,5][which(paramss[,1]=="A3")]), probs = c(0.025,0.975))),nrow = 1),2)
  A3a.expl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,6][which(paramss[,1]=="A3")]), probs = c(0.025,0.975))),nrow = 1),2)
  A3b.expl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,7][which(paramss[,1]=="A3")]), probs = c(0.025,0.975))),nrow = 1),2)
  A3c.expl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,8][which(paramss[,1]=="A3")]), probs = c(0.025,0.975))),nrow = 1),2)

  A4a.ex <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,2][which(paramss[,1]=="A4")]), probs = c(0.025,0.975))),nrow = 1),2)
  A4b.ex <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,3][which(paramss[,1]=="A4")]), probs = c(0.025,0.975))),nrow = 1),2)
  A4a.pl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,4][which(paramss[,1]=="A4")]), probs = c(0.025,0.975))),nrow = 1),2)
  A4b.pl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,5][which(paramss[,1]=="A4")]), probs = c(0.025,0.975))),nrow = 1),2)
  A4a.expl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,6][which(paramss[,1]=="A4")]), probs = c(0.025,0.975))),nrow = 1),2)
  A4b.expl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,7][which(paramss[,1]=="A4")]), probs = c(0.025,0.975))),nrow = 1),2)
  A4c.expl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,8][which(paramss[,1]=="A4")]), probs = c(0.025,0.975))),nrow = 1),2)

  A5a.ex <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,2][which(paramss[,1]=="A5")]), probs = c(0.025,0.975))),nrow = 1),2)
  A5b.ex <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,3][which(paramss[,1]=="A5")]), probs = c(0.025,0.975))),nrow = 1),2)
  A5a.pl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,4][which(paramss[,1]=="A5")]), probs = c(0.025,0.975))),nrow = 1),2)
  A5b.pl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,5][which(paramss[,1]=="A5")]), probs = c(0.025,0.975))),nrow = 1),2)
  A5a.expl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,6][which(paramss[,1]=="A5")]), probs = c(0.025,0.975))),nrow = 1),2)
  A5b.expl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,7][which(paramss[,1]=="A5")]), probs = c(0.025,0.975))),nrow = 1),2)
  A5c.expl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,8][which(paramss[,1]=="A5")]), probs = c(0.025,0.975))),nrow = 1),2)

  A6a.ex <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,2][which(paramss[,1]=="A6")]), probs = c(0.025,0.975))),nrow = 1),2)
  A6b.ex <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,3][which(paramss[,1]=="A6")]), probs = c(0.025,0.975))),nrow = 1),2)
  A6a.pl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,4][which(paramss[,1]=="A6")]), probs = c(0.025,0.975))),nrow = 1),2)
  A6b.pl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,5][which(paramss[,1]=="A6")]), probs = c(0.025,0.975))),nrow = 1),2)
  A6a.expl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,6][which(paramss[,1]=="A6")]), probs = c(0.025,0.975))),nrow = 1),2)
  A6b.expl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,7][which(paramss[,1]=="A6")]), probs = c(0.025,0.975))),nrow = 1),2)
  A6c.expl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,8][which(paramss[,1]=="A6")]), probs = c(0.025,0.975))),nrow = 1),2)

  A7a.ex <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,2][which(paramss[,1]=="A7")]), probs = c(0.025,0.975))),nrow = 1),2)
  A7b.ex <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,3][which(paramss[,1]=="A7")]), probs = c(0.025,0.975))),nrow = 1),2)
  A7a.pl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,4][which(paramss[,1]=="A7")]), probs = c(0.025,0.975))),nrow = 1),2)
  A7b.pl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,5][which(paramss[,1]=="A7")]), probs = c(0.025,0.975))),nrow = 1),2)
  A7a.expl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,6][which(paramss[,1]=="A7")]), probs = c(0.025,0.975))),nrow = 1),2)
  A7b.expl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,7][which(paramss[,1]=="A7")]), probs = c(0.025,0.975))),nrow = 1),2)
  A7c.expl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,8][which(paramss[,1]=="A7")]), probs = c(0.025,0.975))),nrow = 1),2)

  A8a.ex <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,2][which(paramss[,1]=="A8")]), probs = c(0.025,0.975))),nrow = 1),2)
  A8b.ex <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,3][which(paramss[,1]=="A8")]), probs = c(0.025,0.975))),nrow = 1),2)
  A8a.pl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,4][which(paramss[,1]=="A8")]), probs = c(0.025,0.975))),nrow = 1),2)
  A8b.pl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,5][which(paramss[,1]=="A8")]), probs = c(0.025,0.975))),nrow = 1),2)
  A8a.expl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,6][which(paramss[,1]=="A8")]), probs = c(0.025,0.975))),nrow = 1),2)
  A8b.expl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,7][which(paramss[,1]=="A8")]), probs = c(0.025,0.975))),nrow = 1),2)
  A8c.expl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,8][which(paramss[,1]=="A8")]), probs = c(0.025,0.975))),nrow = 1),2)

  A9a.ex <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,2][which(paramss[,1]=="A9")]), probs = c(0.025,0.975))),nrow = 1),2)
  A9b.ex <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,3][which(paramss[,1]=="A9")]), probs = c(0.025,0.975))),nrow = 1),2)
  A9a.pl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,4][which(paramss[,1]=="A9")]), probs = c(0.025,0.975))),nrow = 1),2)
  A9b.pl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,5][which(paramss[,1]=="A9")]), probs = c(0.025,0.975))),nrow = 1),2)
  A9a.expl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,6][which(paramss[,1]=="A9")]), probs = c(0.025,0.975))),nrow = 1),2)
  A9b.expl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,7][which(paramss[,1]=="A9")]), probs = c(0.025,0.975))),nrow = 1),2)
  A9c.expl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,8][which(paramss[,1]=="A9")]), probs = c(0.025,0.975))),nrow = 1),2)

  naa.ex <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,2][which(paramss[,1]=="NA")]), probs = c(0.025,0.975))),nrow = 1),2)
  nab.ex <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,3][which(paramss[,1]=="NA")]), probs = c(0.025,0.975))),nrow = 1),2)
  naa.pl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,4][which(paramss[,1]=="NA")]), probs = c(0.025,0.975))),nrow = 1),2)
  nab.pl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,5][which(paramss[,1]=="NA")]), probs = c(0.025,0.975))),nrow = 1),2)
  naa.expl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,6][which(paramss[,1]=="NA")]), probs = c(0.025,0.975))),nrow = 1),2)
  nab.expl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,7][which(paramss[,1]=="NA")]), probs = c(0.025,0.975))),nrow = 1),2)
  nac.expl <- round(matrix(as.numeric(stats::quantile(
    as.numeric(paramss[,8][which(paramss[,1]=="NA")]), probs = c(0.025,0.975))),nrow = 1),2)

  arche <- noquote(matrix(c("A1","A2","A3","A4","A5","A6","A7","A8","A9",
                            "NA"),ncol=1))

  ##### 95% C.I of parameter estimates for communities under each archetype

  a1 <- matrix(c(A1a.ex[,1],A1a.ex[,2],A1b.ex[,1],A1b.ex[,2],A1a.pl[,1],
                 A1a.pl[,2],A1b.pl[,1],A1b.pl[,2],A1a.expl[,1],
                 A1a.expl[,2],A1b.expl[,1],A1b.expl[,2],A1c.expl[,1],
                 A1c.expl[,2]),ncol=14)

  a2 <- matrix(c(A2a.ex[,1],A2a.ex[,2],A2b.ex[,1],A2b.ex[,2],A2a.pl[,1],
                 A2a.pl[,2],A2b.pl[,1],A2b.pl[,2],A2a.expl[,1],
                 A2a.expl[,2],A2b.expl[,1],A2b.expl[,2],A2c.expl[,1],
                 A2c.expl[,2]),ncol = 14)

  a3 <- matrix(c(A3a.ex[,1],A3a.ex[,2],A3b.ex[,1],A3b.ex[,2],A3a.pl[,1],
                 A3a.pl[,2],A3b.pl[,1],A3b.pl[,2],A3a.expl[,1],
                 A3a.expl[,2],A3b.expl[,1],A3b.expl[,2],A3c.expl[,1],
                 A3c.expl[,2]),ncol = 14)

  a4 <- matrix(c(A4a.ex[,1],A4a.ex[,2],A4b.ex[,1],A4b.ex[,2],A4a.pl[,1],
                 A4a.pl[,2],A4b.pl[,1],A4b.pl[,2],A4a.expl[,1],
                 A4a.expl[,2],A4b.expl[,1],A4b.expl[,2],A4c.expl[,1],
                 A4c.expl[,2]),ncol = 14)

  a5 <- matrix(c(A5a.ex[,1],A5a.ex[,2],A5b.ex[,1],A5b.ex[,2],A5a.pl[,1],
                 A5a.pl[,2],A5b.pl[,1],A5b.pl[,2],A5a.expl[,1],
                 A5a.expl[,2],A5b.expl[,1],A5b.expl[,2],A5c.expl[,1],
                 A5c.expl[,2]),ncol = 14)

  a6 <- matrix(c(A6a.ex[,1],A6a.ex[,2],A6b.ex[,1],A6b.ex[,2],A6a.pl[,1],
                 A6a.pl[,2],A6b.pl[,1],A6b.pl[,2],A6a.expl[,1],
                 A6a.expl[,2],A6b.expl[,1],A6b.expl[,2],A6c.expl[,1],
                 A6c.expl[,2]),ncol = 14)

  a7 <- matrix(c(A7a.ex[,1],A7a.ex[,2],A7b.ex[,1],A7b.ex[,2],A7a.pl[,1],
                 A7a.pl[,2],A7b.pl[,1],A7b.pl[,2],A7a.expl[,1],
                 A7a.expl[,2],A7b.expl[,1],A7b.expl[,2],A7c.expl[,1],
                 A7c.expl[,2]),ncol = 14)

  a8 <- matrix(c(A8a.ex[,1],A8a.ex[,2],A8b.ex[,1],A8b.ex[,2],A8a.pl[,1],
                 A8a.pl[,2],A8b.pl[,1],A8b.pl[,2],A8a.expl[,1],
                 A8a.expl[,2],A8b.expl[,1],A8b.expl[,2],A8c.expl[,1],
                 A8c.expl[,2]),ncol = 14)

  a9 <- matrix(c(A9a.ex[,1],A9a.ex[,2],A9b.ex[,1],A9b.ex[,2],A9a.pl[,1],
                 A9a.pl[,2],A9b.pl[,1],A9b.pl[,2],A9a.expl[,1],
                 A9a.expl[,2],A9b.expl[,1],A9b.expl[,2],A9c.expl[,1],
                 A9c.expl[,2]),ncol = 14)

  na2 <- matrix(c(naa.ex[,1],naa.ex[,2],nab.ex[,1],nab.ex[,2],naa.pl[,1],
                  naa.pl[,2],nab.pl[,1],nab.pl[,2],naa.expl[,1],
                  naa.expl[,2],nab.expl[,1],nab.expl[,2],nac.expl[,1],
                  nac.expl[,2]),ncol = 14)

  N <- matrix(c(length(which(paramss[,1]=="A1")),length(which(paramss[,1]=="A2")),
                length(which(paramss[,1]=="A3")),length(which(paramss[,1]=="A4")),
                length(which(paramss[,1]=="A5")),length(which(paramss[,1]=="A6")),
                length(which(paramss[,1]=="A7")),length(which(paramss[,1]=="A8")),
                length(which(paramss[,1]=="A9")),length(which(paramss[,1]=="NA"))),
              nrow = 10,ncol = 1)

  PARAMSS <- rbind(a1,a2,a3,a4,a5,a6,a7,a8,a9,na2)
  PARAMS <- noquote(cbind(arche[,1],N,PARAMSS))
  paras_c.i <- data.frame(matrix(PARAMS, nrow = nrow(PARAMS), ncol(PARAMS)));paras_c.i

  names(paras_c.i) <- c("Arch","n","L.a_ex","U.a_ex","L.b_ex","U.b_ex",
                        "L.a_pl","U.a_pl","L.b_pl","U.b_pl","L.a_expl",
                        "U.a_expl","L.b_expl","U.b_expl","L.c_expl","U.c_expl")

  #################################################################################
  ############ Regression model comparisons for each Archetype ####################


  ###################################  A1  ########################################

  if(length(which(paramss[,1]=="A1"))==0){
    ex.a1_percent <- NA
    pl.a1_percent <- NA
    ex.pl.a1_percent <- NA

  }else if(length(which(paramss[,1]=="A1"))!=0){
    AICs.a1 <- all.AICs[which(paramss[,1]=="A1")]
    ex.a1 <- c()
    pl.a1 <- c()
    ex.pl.a1 <- c()
    for (jim in 1:length(AICs.a1)) {
      if(AICs.a1[[jim]][,2][1]<AICs.a1[[jim]][,2][2]){
        ex.a1[jim] <- "ex.a1"
      }else if(AICs.a1[[jim]][,2][2]<AICs.a1[[jim]][,2][1]){
        pl.a1[jim] <- "pl.a1"
      }
    }
    for (djm in 1:length(AICs.a1)) {
      if(min(AICs.a1[[djm]][,2][1],AICs.a1[[djm]][,2][2],
             AICs.a1[[djm]][,2][3])==AICs.a1[[djm]][,2][3]){
        ex.pl.a1[djm] <- "ex.pl.a1"
      }
    }

    ex.a1_percent <- ((length(which(ex.a1=="ex.a1")))/length(which(paramss[,1]=="A1")))*100
    pl.a1_percent <- ((length(which(pl.a1=="pl.a1")))/length(which(paramss[,1]=="A1")))*100
    ex.pl.a1_percent <- ((length(which(ex.pl.a1=="ex.pl.a1")))/length(which(
      paramss[,1]=="A1")))*100
  }

  #########################################  A2  #########################################

  if(length(which(paramss[,1]=="A2"))==0){
    ex.a2_percent <- NA
    pl.a2_percent <- NA
    ex.pl.a2_percent <- NA

  }else if(length(which(paramss[,1]=="A2"))!=0){
    AICs.a2 <- all.AICs[which(paramss[,1]=="A2")]
    ex.a2 <- c()
    pl.a2 <- c()
    ex.pl.a2 <- c()
    for (cjs in 1:length(AICs.a2)) {
      if(AICs.a2[[cjs]][,2][1]<AICs.a2[[cjs]][,2][2]){
        ex.a2[cjs] <- "ex.a2"
      }else if(AICs.a2[[cjs]][,2][2]<AICs.a2[[cjs]][,2][1]){
        pl.a2[cjs] <- "pl.a2"
      }
    }
    for (jsc in 1:length(AICs.a2)) {
      if(min(AICs.a2[[jsc]][,2][1],AICs.a2[[jsc]][,2][2],
             AICs.a2[[jsc]][,2][3])==AICs.a2[[jsc]][,2][3]){
        ex.pl.a2[jsc] <- "ex.pl.a2"
      }
    }

    ex.a2_percent <- ((length(which(ex.a2=="ex.a2")))/length(which(paramss[,1]=="A2")))*100
    pl.a2_percent <- ((length(which(pl.a2=="pl.a2")))/length(which(paramss[,1]=="A2")))*100
    ex.pl.a2_percent <- ((length(which(ex.pl.a2=="ex.pl.a2")))/length(which(
      paramss[,1]=="A2")))*100
  }

  ##################################  A3  #########################################

  if(length(which(paramss[,1]=="A3"))==0){
    ex.a3_percent <- NA
    pl.a3_percent <- NA
    ex.pl.a3_percent <- NA

  }else if(length(which(paramss[,1]=="A3"))!=0){
    AICs.a3 <- all.AICs[which(paramss[,1]=="A3")]
    ex.a3 <- c()
    pl.a3 <- c()
    ex.pl.a3 <- c()
    for (cej in 1:length(AICs.a3)) {
      if(AICs.a3[[cej]][,2][1]<AICs.a3[[cej]][,2][2]){
        ex.a3[cej] <- "ex.a3"
      }else if(AICs.a3[[cej]][,2][2]<AICs.a3[[cej]][,2][1]){
        pl.a3[cej] <- "pl.a3"
      }
    }
    for (lik in 1:length(AICs.a3)) {
      if(min(AICs.a3[[lik]][,2][1],AICs.a3[[lik]][,2][2],
             AICs.a3[[lik]][,2][3])==AICs.a3[[lik]][,2][3]){
        ex.pl.a3[lik] <- "ex.pl.a3"
      }
    }

    ex.a3_percent <- ((length(which(ex.a3=="ex.a3")))/length(which(paramss[,1]=="A3")))*100
    pl.a3_percent <- ((length(which(pl.a3=="pl.a3")))/length(which(paramss[,1]=="A3")))*100
    ex.pl.a3_percent <- ((length(which(ex.pl.a3=="ex.pl.a3")))/length(which(
      paramss[,1]=="A3")))*100
  }

  ##################################  A4  #########################################

  if(length(which(paramss[,1]=="A4"))==0){
    ex.a4_percent <- NA
    pl.a4_percent <- NA
    ex.pl.a4_percent <- NA

  }else if(length(which(paramss[,1]=="A4"))!=0){
    AICs.a4 <- all.AICs[which(paramss[,1]=="A4")]
    ex.a4 <- c()
    pl.a4 <- c()
    ex.pl.a4 <- c()
    for (jec in 1:length(AICs.a4)) {
      if(AICs.a4[[jec]][,2][1]<AICs.a4[[jec]][,2][2]){
        ex.a4[jec] <- "ex.a4"
      }else if(AICs.a4[[jec]][,2][2]<AICs.a4[[jec]][,2][1]){
        pl.a4[jec] <- "pl.a4"
      }
    }
    for (kil in 1:length(AICs.a4)) {
      if(min(AICs.a4[[kil]][,2][1],AICs.a4[[kil]][,2][2],
             AICs.a4[[kil]][,2][3])==AICs.a4[[kil]][,2][3]){
        ex.pl.a4[kil] <- "ex.pl.a4"
      }
    }

    ex.a4_percent <- ((length(which(ex.a4=="ex.a4")))/length(which(paramss[,1]=="A4")))*100
    pl.a4_percent <- ((length(which(pl.a4=="pl.a4")))/length(which(paramss[,1]=="A4")))*100
    ex.pl.a4_percent <- ((length(which(ex.pl.a4=="ex.pl.a4")))/length(which(
      paramss[,1]=="A4")))*100
  }

  #########################################  A5  #########################################

  if(length(which(paramss[,1]=="A5"))==0){
    ex.a5_percent <- NA
    pl.a5_percent <- NA
    ex.pl.a5_percent <- NA

  }else if(length(which(paramss[,1]=="A5"))!=0){
    AICs.a5 <- all.AICs[which(paramss[,1]=="A5")]
    ex.a5 <- c()
    pl.a5 <- c()
    ex.pl.a5 <- c()
    for (jel in 1:length(AICs.a5)) {
      if(AICs.a5[[jel]][,2][1]<AICs.a5[[jel]][,2][2]){
        ex.a5[jel] <- "ex.a5"
      }else if(AICs.a5[[jel]][,2][2]<AICs.a5[[jel]][,2][1]){
        pl.a5[jel] <- "pl.a5"
      }
    }
    for (kic in 1:length(AICs.a5)) {
      if(min(AICs.a5[[kic]][,2][1],AICs.a5[[kic]][,2][2],
             AICs.a5[[kic]][,2][3])==AICs.a5[[kic]][,2][3]){
        ex.pl.a5[kic] <- "ex.pl.a5"
      }
    }

    ex.a5_percent <- ((length(which(ex.a5=="ex.a5")))/length(which(paramss[,1]=="A5")))*100
    pl.a5_percent <- ((length(which(pl.a5=="pl.a5")))/length(which(paramss[,1]=="A5")))*100
    ex.pl.a5_percent <- ((length(which(ex.pl.a5=="ex.pl.a5")))/length(
      which(paramss[,1]=="A5")))*100
  }

  #########################################  A6  #########################################

  if(length(which(paramss[,1]=="A6"))==0){
    ex.a6_percent <- NA
    pl.a6_percent <- NA
    ex.pl.a6_percent <- NA

  }else if(length(which(paramss[,1]=="A6"))!=0){
    AICs.a6 <- all.AICs[which(paramss[,1]=="A6")]
    ex.a6 <- c()
    pl.a6 <- c()
    ex.pl.a6 <- c()
    for (jc in 1:length(AICs.a6)) {
      if(AICs.a6[[jc]][,2][1]<AICs.a6[[jc]][,2][2]){
        ex.a6[jc] <- "ex.a6"
      }else if(AICs.a6[[jc]][,2][2]<AICs.a6[[jc]][,2][1]){
        pl.a6[jc] <- "pl.a6"
      }
    }
    for (kc in 1:length(AICs.a6)) {
      if(min(AICs.a6[[kc]][,2][1],AICs.a6[[kc]][,2][2],
             AICs.a6[[kc]][,2][3])==AICs.a6[[kc]][,2][3]){
        ex.pl.a6[kc] <- "ex.pl.a6"
      }
    }

    ex.a6_percent <- ((length(which(ex.a6=="ex.a6")))/length(which(paramss[,1]=="A6")))*100
    pl.a6_percent <- ((length(which(pl.a6=="pl.a6")))/length(which(paramss[,1]=="A6")))*100
    ex.pl.a6_percent <- ((length(which(ex.pl.a6=="ex.pl.a6")))/length(which(
      paramss[,1]=="A6")))*100
  }

  #########################################  A7  #########################################

  if(length(which(paramss[,1]=="A7"))==0){
    ex.a7_percent <- NA
    pl.a7_percent <- NA
    ex.pl.a7_percent <- NA

  }else if(length(which(paramss[,1]=="A7"))!=0){
    AICs.a7 <- all.AICs[which(paramss[,1]=="A7")]
    ex.a7 <- c()
    pl.a7 <- c()
    ex.pl.a7 <- c()
    for (jl in 1:length(AICs.a7)) {
      if(AICs.a7[[jl]][,2][1]<AICs.a7[[jl]][,2][2]){
        ex.a7[jl] <- "ex.a7"
      }else if(AICs.a7[[jl]][,2][2]<AICs.a7[[jl]][,2][1]){
        pl.a7[jl] <- "pl.a7"
      }
    }
    for (kl in 1:length(AICs.a7)) {
      if(min(AICs.a7[[kl]][,2][1],AICs.a7[[kl]][,2][2],
             AICs.a7[[kl]][,2][3])==AICs.a7[[kl]][,2][3]){
        ex.pl.a7[kl] <- "ex.pl.a7"
      }
    }

    ex.a7_percent <- ((length(which(ex.a7=="ex.a7")))/length(which(paramss[,1]=="A7")))*100
    pl.a7_percent <- ((length(which(pl.a7=="pl.a7")))/length(which(paramss[,1]=="A7")))*100
    ex.pl.a7_percent <- ((length(which(ex.pl.a7=="ex.pl.a7")))/length(which(
      paramss[,1]=="A7")))*100
  }

  ##########################################  A8  #########################################

  if(length(which(paramss[,1]=="A8"))==0){
    ex.a8_percent <- NA
    pl.a8_percent <- NA
    ex.pl.a8_percent <- NA

  }else if(length(which(paramss[,1]=="A8"))!=0){
    AICs.a8 <- all.AICs[which(paramss[,1]=="A8")]
    ex.a8 <- c()
    pl.a8 <- c()
    ex.pl.a8 <- c()
    for (jv in 1:length(AICs.a8)) {
      if(AICs.a8[[jv]][,2][1]<AICs.a8[[jv]][,2][2]){
        ex.a8[jv] <- "ex.a8"
      }else if(AICs.a8[[jv]][,2][2]<AICs.a8[[jv]][,2][1]){
        pl.a8[jv] <- "pl.a8"
      }
    }
    for (kv in 1:length(AICs.a8)) {
      if(min(AICs.a8[[kv]][,2][1],AICs.a8[[kv]][,2][2],
             AICs.a8[[kv]][,2][3])==AICs.a8[[kv]][,2][3]){
        ex.pl.a8[kv] <- "ex.pl.a8"
      }
    }

    ex.a8_percent <- ((length(which(ex.a8=="ex.a8")))/length(which(paramss[,1]=="A8")))*100
    pl.a8_percent <- ((length(which(pl.a8=="pl.a8")))/length(which(paramss[,1]=="A8")))*100
    ex.pl.a8_percent <- ((length(which(ex.pl.a8=="ex.pl.a8")))/length(which(
      paramss[,1]=="A8")))*100
  }

  #########################################  A9  #########################################

  if(length(which(paramss[,1]=="A9"))==0){
    ex.a9_percent <- NA
    pl.a9_percent <- NA
    ex.pl.a9_percent <- NA

  }else if(length(which(paramss[,1]=="A9"))!=0){
    AICs.a9 <- all.AICs[which(paramss[,1]=="A9")]
    ex.a9 <- c()
    pl.a9 <- c()
    ex.pl.a9 <- c()
    for (jj in 1:length(AICs.a9)) {
      if(AICs.a9[[jj]][,2][1]<AICs.a9[[jj]][,2][2]){
        ex.a9[jj] <- "ex.a9"
      }else if(AICs.a9[[jj]][,2][2]<AICs.a9[[jj]][,2][1]){
        pl.a9[jj] <- "pl.a9"
      }
    }
    for (kk in 1:length(AICs.a9)) {
      if(min(AICs.a9[[kk]][,2][1],AICs.a9[[kk]][,2][2],
             AICs.a9[[kk]][,2][3])==AICs.a9[[kk]][,2][3]){
        ex.pl.a9[kk] <- "ex.pl.a9"
      }
    }

    ex.a9_percent <- ((length(which(ex.a9=="ex.a9")))/length(which(paramss[,1]=="A9")))*100
    pl.a9_percent <- ((length(which(pl.a9=="pl.a9")))/length(which(paramss[,1]=="A9")))*100
    ex.pl.a9_percent <- ((length(which(ex.pl.a9=="ex.pl.a9")))/length(which(
      paramss[,1]=="A9")))*100
  }

  ############################################  NA  ######################################

  if(length(which(paramss[,1]=="NA"))==0){
    ex.na_percent <- NA
    pl.na_percent <- NA
    ex.pl.na_percent <- NA

  }else if(length(which(paramss[,1]=="NA"))!=0){
    AICs.na <- all.AICs[which(paramss[,1]=="NA")]
    ex.na <- c()
    pl.na <- c()
    ex.pl.na <- c()
    for (ji in 1:length(AICs.na)) {
      if(AICs.na[[ji]][,2][1]<AICs.na[[ji]][,2][2]){
        ex.na[ji] <- "ex.na"
      }else if(AICs.na[[ji]][,2][2]<AICs.na[[ji]][,2][1]){
        pl.na[ji] <- "pl.na"
      }
    }
    for (ik in 1:length(AICs.na)) {
      if(min(AICs.na[[ik]][,2][1],AICs.na[[ik]][,2][2],
             AICs.na[[ik]][,2][3])==AICs.na[[ik]][,2][3]){
        ex.pl.na[ik] <- "ex.pl.na"
      }
    }

    ex.na_percent <- ((length(which(ex.na=="ex.na")))/length(which(paramss[,1]=="NA")))*100
    pl.na_percent <- ((length(which(pl.na=="pl.na")))/length(which(paramss[,1]=="NA")))*100
    ex.pl.na_percent <- ((length(which(ex.pl.na=="ex.pl.na")))/length(
      which(paramss[,1]=="NA")))*100

  }
  ########################################################################################
  ############################# Params_c.i with archetype % ##############################

  ex.perc <- round(matrix(c(ex.a1_percent, ex.a2_percent, ex.a3_percent,
                            ex.a4_percent, ex.a5_percent, ex.a6_percent,
                            ex.a7_percent, ex.a8_percent, ex.a9_percent,
                            ex.na_percent), nrow = length(paras_c.i[,1]),
                          ncol = 1),1)

  ex.perc[,1][is.nan(ex.perc[,1])] <- NA

  pl.perc <- round(matrix(c(pl.a1_percent, pl.a2_percent, pl.a3_percent,
                            pl.a4_percent, pl.a5_percent, pl.a6_percent,
                            pl.a7_percent, pl.a8_percent, pl.a9_percent,
                            pl.na_percent), nrow = length(paras_c.i[,1]),
                          ncol = 1),1)

  pl.perc[,1][is.nan(pl.perc[,1])] <- NA

  ex.pl.perc <- round(matrix(c(ex.pl.a1_percent, ex.pl.a2_percent,
                               ex.pl.a3_percent, ex.pl.a4_percent,
                               ex.pl.a5_percent, ex.pl.a6_percent,
                               ex.pl.a7_percent, ex.pl.a8_percent,
                               ex.pl.a9_percent, ex.pl.na_percent),
                             nrow = length(paras_c.i[,1]), ncol = 1),1)

  ex.pl.perc[,1][is.nan(ex.pl.perc[,1])] <- NA;

  paras_c.i <- tibble::add_column(paras_c.i, "ex_%" = ex.perc[,1], .after="n")
  paras_c.i <- tibble::add_column(paras_c.i,"p.l_%" = pl.perc[,1],
                                  .after="U.b_ex")
  paras_c.i <- tibble::add_column(paras_c.i,"ex.pl_%" = ex.pl.perc[,1],
                                  .after="U.b_pl")

  paras_c.i <- tidyr::unite(paras_c.i, col = "a.ex", c("L.a_ex","U.a_ex"), sep = ",")
  paras_c.i <- tidyr::unite(paras_c.i, col = "b.ex", c("L.b_ex","U.b_ex"), sep = ",")
  paras_c.i <- tidyr::unite(paras_c.i, col = "a.pl", c("L.a_pl","U.a_pl"), sep = ",")
  paras_c.i <- tidyr::unite(paras_c.i, col = "b.pl", c("L.b_pl","U.b_pl"), sep = ",")
  paras_c.i <- tidyr::unite(paras_c.i, col = "a.expl", c("L.a_expl","U.a_expl"),
                            sep = ",")
  paras_c.i <- tidyr::unite(paras_c.i, col = "b.expl", c("L.b_expl","U.b_expl"),
                            sep = ",")
  paras_c.i <- tidyr::unite(paras_c.i, col = "c.expl", c("L.c_expl","U.c_expl"),
                            sep = ",")

  paras_c.i <- paras_c.i[stats::complete.cases(paras_c.i),]

  for(vcj in (4:ncol(paras_c.i))[c(-3,-6)]){
    paras_c.i[,vcj] <- sapply(1:nrow(paras_c.i), function(z) do.call(
      sprintf, as.list(c("[%s%s", c(paras_c.i[,vcj][z], "]")))))
  }

  ########################################################################################
  ####################################### rsq Summary #####################################
  ########################################################################################

  ######################################### All.rsq #######################################
  ex.r2 <- c()
  pl.r2 <- c()
  ex.pl.r2 <- c()
  for (i in 1:length(r2)) {
    if(r2[[i]][1]>0.95){
      ex.r2[i] <- "ex.r2"
    }
  }
  for (i in 1:length(r2)) {
    if(r2[[i]][2]>0.95){
      pl.r2[i] <- "pl.r2"
    }
  }
  for (j in 1:length(r2)) {
    if(r2[[j]][3]>0.95){
      ex.pl.r2[j] <- "ex.pl.r2"
    }
  }

  ###################################### A1 ##################################

  if(length(which(Archs=="A1"))==0){

    ex.r2.a1 <- NA
    pl.r2.a1 <- NA
    ex.pl.r2.a1 <- NA

  }else if(length(which(Archs=="A1"))!=0){
    r2.a1 <- r2[which(Archs=="A1")]
    ex.r2.a1 <- c()
    pl.r2.a1 <- c()
    ex.pl.r2.a1 <- c()
    for (i in 1:length(r2.a1)) {
      if(r2.a1[[i]][1]>0.95){
        ex.r2.a1[i] <- "ex.r2.a1"
      }
    }
    for (i in 1:length(r2.a1)) {
      if(r2.a1[[i]][2]>0.95){
        pl.r2.a1[i] <- "pl.r2.a1"
      }
    }
    for (j in 1:length(r2.a1)) {
      if(r2.a1[[j]][3]>0.95){
        ex.pl.r2.a1[j] <- "ex.pl.r2.a1"
      }
    }
  }

  #######################################  A2  ####################################

  if(length(which(Archs=="A2"))==0){

    ex.r2.a2 <- NA
    pl.r2.a2 <- NA
    ex.pl.r2.a2 <- NA

  }else if(length(which(Archs=="A2"))!=0){
    r2.a2 <- r2[which(Archs=="A2")]
    ex.r2.a2 <- c()
    pl.r2.a2 <- c()
    ex.pl.r2.a2 <- c()
    for (i in 1:length(r2.a2)) {
      if(r2.a2[[i]][1]>0.95){
        ex.r2.a2[i] <- "ex.r2.a2"
      }
    }
    for (i in 1:length(r2.a2)) {
      if(r2.a2[[i]][2]>0.95){
        pl.r2.a2[i] <- "pl.r2.a2"
      }
    }
    for (j in 1:length(r2.a2)) {
      if(r2.a2[[j]][3]>0.95){
        ex.pl.r2.a2[j] <- "ex.pl.r2.a2"
      }
    }
  }

  ######################################  A3  #####################################

  if(length(which(Archs=="A3"))==0){

    ex.r2.a3 <- NA
    pl.r2.a3 <- NA
    ex.pl.r2.a3 <- NA

  }else if(length(which(Archs=="A3"))!=0){
    r2.a3 <- r2[which(Archs=="A3")]
    ex.r2.a3 <- c()
    pl.r2.a3 <- c()
    ex.pl.r2.a3 <- c()
    for (i in 1:length(r2.a3)) {
      if(r2.a3[[i]][1]>0.95){
        ex.r2.a3[i] <- "ex.r2.a3"
      }
    }
    for (i in 1:length(r2.a3)) {
      if(r2.a3[[i]][2]>0.95){
        pl.r2.a3[i] <- "pl.r2.a3"
      }
    }
    for (j in 1:length(r2.a3)) {
      if(r2.a3[[j]][3]>0.95){
        ex.pl.r2.a3[j] <- "ex.pl.r2.a3"
      }
    }
  }

  ######################################  A4  #####################################

  if(length(which(Archs=="A4"))==0){

    ex.r2.a4 <- NA
    pl.r2.a4 <- NA
    ex.pl.r2.a4 <- NA

  }else if(length(which(Archs=="A4"))!=0){
    r2.a4 <- r2[which(Archs=="A4")]
    ex.r2.a4 <- c()
    pl.r2.a4 <- c()
    ex.pl.r2.a4 <- c()
    for (i in 1:length(r2.a4)) {
      if(r2.a4[[i]][1]>0.95){
        ex.r2.a4[i] <- "ex.r2.a4"
      }
    }
    for (i in 1:length(r2.a4)) {
      if(r2.a4[[i]][2]>0.95){
        pl.r2.a4[i] <- "pl.r2.a4"
      }
    }
    for (j in 1:length(r2.a4)) {
      if(r2.a4[[j]][3]>0.95){
        ex.pl.r2.a4[j] <- "ex.pl.r2.a4"
      }
    }
  }

  ######################################  A5  #####################################

  if(length(which(Archs=="A5"))==0){

    ex.r2.a5 <- NA
    pl.r2.a5 <- NA
    ex.pl.r2.a5 <- NA

  }else if(length(which(Archs=="A5"))!=0){
    r2.a5 <- r2[which(Archs=="A5")]
    ex.r2.a5 <- c()
    pl.r2.a5 <- c()
    ex.pl.r2.a5 <- c()
    for (i in 1:length(r2.a5)) {
      if(r2.a5[[i]][1]>0.95){
        ex.r2.a5[i] <- "ex.r2.a5"
      }
    }
    for (i in 1:length(r2.a5)) {
      if(r2.a5[[i]][2]>0.95){
        pl.r2.a5[i] <- "pl.r2.a5"
      }
    }
    for (j in 1:length(r2.a5)) {
      if(r2.a5[[j]][3]>0.95){
        ex.pl.r2.a5[j] <- "ex.pl.r2.a5"
      }
    }
  }

  ######################################  A6  #####################################

  if(length(which(Archs=="A6"))==0){

    ex.r2.a6 <- NA
    pl.r2.a6 <- NA
    ex.pl.r2.a6 <- NA

  }else if(length(which(Archs=="A6"))!=0){
    r2.a6 <- r2[which(Archs=="A6")]
    ex.r2.a6 <- c()
    pl.r2.a6 <- c()
    ex.pl.r2.a6 <- c()
    for (i in 1:length(r2.a6)) {
      if(r2.a6[[i]][1]>0.95){
        ex.r2.a6[i] <- "ex.r2.a6"
      }
    }
    for (i in 1:length(r2.a6)) {
      if(r2.a6[[i]][2]>0.95){
        pl.r2.a6[i] <- "pl.r2.a6"
      }
    }
    for (j in 1:length(r2.a6)) {
      if(r2.a6[[j]][3]>0.95){
        ex.pl.r2.a6[j] <- "ex.pl.r2.a6"
      }
    }
  }

  ######################################  A7  #####################################

  if(length(which(Archs=="A7"))==0){

    ex.r2.a7 <- NA
    pl.r2.a7 <- NA
    ex.pl.r2.a7 <- NA

  }else if(length(which(Archs=="A7"))!=0){
    r2.a7 <- r2[which(Archs=="A7")]
    ex.r2.a7 <- c()
    pl.r2.a7 <- c()
    ex.pl.r2.a7 <- c()
    for (i in 1:length(r2.a7)) {
      if(r2.a7[[i]][1]>0.95){
        ex.r2.a7[i] <- "ex.r2.a7"
      }
    }
    for (i in 1:length(r2.a7)) {
      if(r2.a7[[i]][2]>0.95){
        pl.r2.a7[i] <- "pl.r2.a7"
      }
    }
    for (j in 1:length(r2.a7)) {
      if(r2.a7[[j]][3]>0.95){
        ex.pl.r2.a7[j] <- "ex.pl.r2.a7"
      }
    }
  }

  ######################################  A8  #####################################

  if(length(which(Archs=="A8"))==0){

    ex.r2.a8 <- NA
    pl.r2.a8 <- NA
    ex.pl.r2.a8 <- NA

  }else if(length(which(Archs=="A8"))!=0){
    r2.a8 <- r2[which(Archs=="A8")]
    ex.r2.a8 <- c()
    pl.r2.a8 <- c()
    ex.pl.r2.a8 <- c()
    for (i in 1:length(r2.a8)) {
      if(r2.a8[[i]][1]>0.95){
        ex.r2.a8[i] <- "ex.r2.a8"
      }
    }
    for (i in 1:length(r2.a8)) {
      if(r2.a8[[i]][2]>0.95){
        pl.r2.a8[i] <- "pl.r2.a8"
      }
    }
    for (j in 1:length(r2.a8)) {
      if(r2.a8[[j]][3]>0.95){
        ex.pl.r2.a8[j] <- "ex.pl.r2.a8"
      }
    }
  }

  ######################################  A9  #####################################

  if(length(which(Archs=="A9"))==0){

    ex.r2.a9 <- NA
    pl.r2.a9 <- NA
    ex.pl.r2.a9 <- NA

  }else if(length(which(Archs=="A9"))!=0){
    r2.a9 <- r2[which(Archs=="A9")]
    ex.r2.a9 <- c()
    pl.r2.a9 <- c()
    ex.pl.r2.a9 <- c()
    for (i in 1:length(r2.a9)) {
      if(r2.a9[[i]][1]>0.95){
        ex.r2.a9[i] <- "ex.r2.a9"
      }
    }
    for (i in 1:length(r2.a9)) {
      if(r2.a9[[i]][2]>0.95){
        pl.r2.a9[i] <- "pl.r2.a9"
      }
    }
    for (j in 1:length(r2.a9)) {
      if(r2.a9[[j]][3]>0.95){
        ex.pl.r2.a9[j] <- "ex.pl.r2.a9"
      }
    }
  }

  ######################################  NA  #####################################

  if(length(which(Archs=="NA"))==0){

    ex.r2.na <- NA
    pl.r2.na <- NA
    ex.pl.r2.na <- NA

  }else if(length(which(Archs=="NA"))!=0){
    r2.na <- r2[which(Archs=="NA")]
    ex.r2.na <- c()
    pl.r2.na <- c()
    ex.pl.r2.na <- c()
    for (i in 1:length(r2.na)) {
      if(r2.na[[i]][1]>0.95){
        ex.r2.na[i] <- "ex.r2.na"
      }
    }
    for (i in 1:length(r2.na)) {
      if(r2.na[[i]][2]>0.95){
        pl.r2.na[i] <- "pl.r2.na"
      }
    }
    for (j in 1:length(r2.na)) {
      if(r2.na[[j]][3]>0.95){
        ex.pl.r2.na[j] <- "ex.pl.r2.na"
      }
    }
  }

  rsq.ex <- matrix(c(length(which(ex.r2.a1=="ex.r2.a1")), length(which(ex.r2.a2=="ex.r2.a2")),
                    length(which(ex.r2.a3=="ex.r2.a3")), length(which(ex.r2.a4=="ex.r2.a4")),
                    length(which(ex.r2.a5=="ex.r2.a5")), length(which(ex.r2.a6=="ex.r2.a6")),
                    length(which(ex.r2.a7=="ex.r2.a7")), length(which(ex.r2.a8=="ex.r2.a8")),
                    length(which(ex.r2.a9=="ex.r2.a9")),length(which(ex.r2.na=="ex.r2.na"))),
                  ncol = 1)
  rsq.pl <- matrix(c(length(which(pl.r2.a1=="pl.r2.a1")), length(which(pl.r2.a2=="pl.r2.a2")),
                    length(which(pl.r2.a3=="pl.r2.a3")), length(which(pl.r2.a4=="pl.r2.a4")),
                    length(which(pl.r2.a5=="pl.r2.a5")), length(which(pl.r2.a6=="pl.r2.a6")),
                    length(which(pl.r2.a7=="pl.r2.a7")), length(which(pl.r2.a8=="pl.r2.a8")),
                    length(which(pl.r2.a9=="pl.r2.a9")),length(which(pl.r2.na=="pl.r2.na"))),
                  ncol = 1)
  rsq.ex.pl <- matrix(c(length(which(ex.pl.r2.a1=="ex.pl.r2.a1")), length(which(
    ex.pl.r2.a2=="ex.pl.r2.a2")),
    length(which(ex.pl.r2.a3=="ex.pl.r2.a3")), length(which(ex.pl.r2.a4=="ex.pl.r2.a4")),
    length(which(ex.pl.r2.a5=="ex.pl.r2.a5")), length(which(
      ex.pl.r2.a6=="ex.pl.r2.a6")),
    length(which(ex.pl.r2.a7=="ex.pl.r2.a7")), length(which(
      ex.pl.r2.a8=="ex.pl.r2.a8")),
    length(which(ex.pl.r2.a9=="ex.pl.r2.a9")),length(which(
      ex.pl.r2.na=="ex.pl.r2.na"))),
    ncol = 1)

  Archse <- noquote(matrix(c("A1","A2","A3","A4","A5","A6","A7","A8","A9",
                            "NA"),ncol=1))
  N <- matrix(c(length(which(Archs=="A1")),length(which(Archs=="A2")),length(which(Archs=="A3")),
                length(which(Archs=="A4")),length(which(Archs=="A5")),length(which(Archs=="A6")),
                length(which(Archs=="A7")),length(which(Archs=="A8")),length(which(Archs=="A9")),
                length(which(Archs=="NA"))),nrow = 10,ncol = 1)


  rsq <- data.frame(Archse[,1], N, rsq.ex, rsq.pl, rsq.ex.pl)
  names(rsq) <- c("Archs","n.a","rsq.ex", "rsq.pl", "rsq.ex.pl")

  all.rsq <- as.data.frame(matrix(c("rsq",length(my.files),length(which(ex.r2=="ex.r2")),
                                   length(which(pl.r2=="pl.r2")),
                                   length(which(ex.pl.r2=="ex.pl.r2"))),
                                 nrow = 1, ncol = 5))

  names(all.rsq) <- c("","n", "Ex", "Pl", "Ex.pl")
  all.rsq <- noquote(unlist(all.rsq))

  r2.s <- list()
  r2.s$rsq.per.Archs <- rsq[which(rowSums(rsq[,2:5]) > 0),]
  r2.s$rsq.all.Communities <- all.rsq

  ####################################### END ######################################################

  myres <-list()
  if(Archetypes == TRUE){
    myres$Archs <- `names<-`(nm_arch, my.files)
  }
  if(AICs == TRUE){
    myres$all.AICs <- `names<-`(all.AICs, my.files)
  }
  if(my.r2 == TRUE){
    myres$r2 <- `names<-`(r2, my.files)
  }
  if(my.r2.s == TRUE){
    myres$r2.s <- r2.s
  }
  if(best.mod2 == TRUE){
    myres$best.mod2 <- both.models
  }
  if(best.mod3 == TRUE){
    myres$best.mod3 <- three.models
  }
  if(params == TRUE){
    myres$params <- paramss
  }
  if(params_c.i == TRUE){
    myres$params_c.i <- paras_c.i
  }
  myres$m.Jo.plots <- print(noquote("Check msco's 'ms' folder in your R version's directory for a 'm.Jo.plots.pdf' file."))

  return(myres)

}


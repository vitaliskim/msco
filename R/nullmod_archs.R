########################################################################################################
#########################                                ###############################################
######################### To replicate Fig. 3 (Jo paper) ###############################################
#########################                                ###############################################
########################################################################################################

nullmod_archs <- function(){
  a1 <- utils::read.csv(system.file("extdata/myCSVs", "133.csv", package = "msco"))
  j.en1 <- Jo.eng(a1, algo="sim2", metric = "raw", nReps = 999,
                  dig = 3, s.dplot = FALSE, All.plots = FALSE, Jo.coeff = FALSE,
                  my.AIC = FALSE, my.rsq = FALSE, Exp_Reg = FALSE, P.law_Reg = FALSE,
                  Exp_p.l_Reg = FALSE, Obs.data = FALSE, Sim.data = FALSE,
                  Jo_val.sim = FALSE, C.I_Jo_val.sim = FALSE, Jo_val.obs = FALSE,
                  Metric = FALSE, Algorithm = FALSE, S.order = FALSE,
                  nmod_stats = FALSE, Pt_Arch_Vals = FALSE, Atype = TRUE,
                  p.n.plot = FALSE, trans = FALSE, lab=FALSE, leg=FALSE, m.n.plot = TRUE)

  a2 <- utils::read.csv(system.file("extdata/myCSVs", "226.csv", package = "msco"))
  j.en2 <- Jo.eng(a2, algo="sim2", metric = "raw", nReps = 999,
                  dig = 3, s.dplot = FALSE, All.plots = FALSE, Jo.coeff = FALSE,
                  my.AIC = FALSE, my.rsq = FALSE, Exp_Reg = FALSE, P.law_Reg = FALSE,
                  Exp_p.l_Reg = FALSE, Obs.data = FALSE, Sim.data = FALSE,
                  Jo_val.sim = FALSE, C.I_Jo_val.sim = FALSE, Jo_val.obs = FALSE,
                  Metric = FALSE, Algorithm = FALSE, S.order = FALSE,
                  nmod_stats = FALSE, Pt_Arch_Vals = FALSE, Atype = TRUE,
                  p.n.plot = FALSE, trans = FALSE, lab=FALSE, leg=FALSE, m.n.plot = TRUE)

  a5 <- utils::read.csv(system.file("extdata/myCSVs", "108.csv", package = "msco"))
  j.en5 <- Jo.eng(a5, algo="sim2", metric = "raw", nReps = 999,
                  dig = 3, s.dplot = FALSE, All.plots = FALSE, Jo.coeff = FALSE,
                  my.AIC = FALSE, my.rsq = FALSE, Exp_Reg = FALSE, P.law_Reg = FALSE,
                  Exp_p.l_Reg = FALSE, Obs.data = FALSE, Sim.data = FALSE,
                  Jo_val.sim = FALSE, C.I_Jo_val.sim = FALSE, Jo_val.obs = FALSE,
                  Metric = FALSE, Algorithm = FALSE, S.order = FALSE,
                  nmod_stats = FALSE, Pt_Arch_Vals = FALSE, Atype = TRUE,
                  p.n.plot = FALSE, trans = FALSE, lab=FALSE, leg=FALSE, m.n.plot = TRUE)

  a9 <- utils::read.csv(system.file("extdata/myCSVs", "274.csv", package = "msco"))
  j.en9 <- Jo.eng(a9, algo="sim2", metric = "raw", nReps = 999,
                  dig = 3, s.dplot = FALSE, All.plots = FALSE, Jo.coeff = FALSE,
                  my.AIC = FALSE, my.rsq = FALSE, Exp_Reg = FALSE, P.law_Reg = FALSE,
                  Exp_p.l_Reg = FALSE, Obs.data = FALSE, Sim.data = FALSE,
                  Jo_val.sim = FALSE, C.I_Jo_val.sim = FALSE, Jo_val.obs = FALSE,
                  Metric = FALSE, Algorithm = FALSE, S.order = FALSE,
                  nmod_stats = FALSE, Pt_Arch_Vals = FALSE, Atype = TRUE,
                  p.n.plot = FALSE, trans = FALSE, lab=FALSE, leg=FALSE, m.n.plot = TRUE)

  a4 <- utils::read.csv(system.file("extdata/myCSVs", "251.csv", package = "msco"))
  j.en4 <- Jo.eng(a4, algo="sim2", metric = "raw", nReps = 999,
                  dig = 3, s.dplot = FALSE, All.plots = FALSE, Jo.coeff = FALSE,
                  my.AIC = FALSE, my.rsq = FALSE, Exp_Reg = FALSE, P.law_Reg = FALSE,
                  Exp_p.l_Reg = FALSE, Obs.data = FALSE, Sim.data = FALSE,
                  Jo_val.sim = FALSE, C.I_Jo_val.sim = FALSE, Jo_val.obs = FALSE,
                  Metric = FALSE, Algorithm = FALSE, S.order = FALSE,
                  nmod_stats = FALSE, Pt_Arch_Vals = FALSE, Atype = TRUE,
                  p.n.plot = FALSE, trans = FALSE, lab=FALSE, leg=FALSE, m.n.plot = TRUE)



  all.plots <- cowplot::ggdraw() +
    cowplot::draw_plot(j.en1$m.n.plot, x = 0.05, y = 0.50, width = .22, height = .5)+
    cowplot::draw_plot(j.en2$m.n.plot, x = 0.25, y = 0.50, width = .22, height = .5)+
    cowplot::draw_plot(j.en5$m.n.plot, x = 0.05, y = 0.05, width = .22, height = .5)+
    cowplot::draw_plot(j.en9$m.n.plot, x = 0.25, y = 0.05, width = .22, height = .5)+
    cowplot::draw_plot(j.en4$m.n.plot, x = 0.45, y = 0.00, width = 0.42, height = 1)+
    cowplot::draw_plot_label(label = c(j.en1$Archetype, j.en2$Archetype, j.en5$Archetype,
                                       j.en9$Archetype, j.en4$Archetype),
                             size = 12, x = c(0.23, 0.43, 0.23, 0.43, 0.83),
                             y = c(0.95, 0.95, 0.5, 0.5, 0.95)) +

    cowplot::draw_image(system.file("logos", "ylab.jpg", package = "msco"),
                        x=0.02, y=0.09, scale = 0.65,
                        width = 0.05) +

    cowplot::draw_image(system.file("logos", "xlab.jpg", package = "msco"),
                        x=0.0, y=0.02, scale = 0.57,
                        height = 0.07)+

    cowplot::draw_image(system.file("logos", "key2.nm.jpg", package = "msco"),
                        x=0.42, y=0.07,scale = 0.15)
  return(all.plots)
}
# st <- Sys.time()
# RNGkind(sample.kind = "Rejection")
# set.seed(14)
# grDevices::pdf(file = paste0(system.file("ms", package = "msco"), "/real.arch.plots.pdf"), paper="a4r", height = 8.27, width = 11.69)
# msco:::nullmod_archs()
# grDevices::dev.off()
# system(paste0('open "', paste0(system.file("ms", package = "msco"), "/real.arch.plots.pdf"), '"'))
# et <- Sys.time();et-st



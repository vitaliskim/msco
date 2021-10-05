#' Archetypes of species co-occurrence patterns
#'
#' This function performs the null model testing of species co-occurrence patterns and generates
#'  the archetypes (e.g. aggregated, segregated, stochastic, etc) of the same patterns using the
#'   species-by-site presence/absence `.csv` data matrices as proxies for ecological communities.
#'    From these archetypes and figure \link[msco]{Arch_schem}, inferences on the forces
#'     structuralising ecological communities can be made.
#'
#'@details We use "archetypes" here to mean different types of patterns of species co-occurrences
#' in communities. These patterns can be clustered, segregated, random, or somewhere in between.
#' See \link[msco]{Arch_schem} function for a schematic diagram of all the possible archetypes
#'  and the discussion on what they mean ecologically.
#'
#' The `Archetypes` function is useful when analyzing multiple species-by-site presence/absence
#'  data matrices at once. If one community matrix is analyzed, the `Archetype` output of the
#'   function \link[msco]{Jo.eng} should suffice.
#'
#' @param my.files A vector containing names of species-by-site presence/absence `.csv` data matrices.
#'  The data matrices should be saved in the working directory.
#' @param algo Simulation algorithm used. The possible options to choose from are: `sim1`,
#'  `sim2`, `sim3`, `sim4`, `sim5`, `sim6`, `sim7`, `sim8`, and `sim9`, all from
#'   Gotelli (2000). `sim2` is highly recommended (see Lagat *et al.,* 202X).
#' @param metric Metric used to quantify the patterns in `s.data`. It has to be multi-species
#'  co-occurrence index (see \link[msco]{j.occ}).
#' @param nReps Number of simulations used in the null model test.
#' @return Returns a `data.frame` with:
#' \itemize{
#' \item names of the species-by-site presence/absence `.csv` data matrices as the `colnames` of the
#'  `data.frame`, and
#' \item the only row of the `data.frame` containing the archetypes of the patterns of species
#'  co-occurrences in ecological communities/matrices (`my.files`). These archetypes must be \eqn{\in \{}"A1",
#'   "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9"\eqn{\}} or "NA". "NA" could be the
#'    combinations of two or more of the nine expected archetypes.
#' }
#' @references
#' \enumerate{
#' \item{Lagat, V. K., Latombe, G. and Hui, C. (202X). *A multi-species co-occurrence index to avoid type
#'  II errors in null model testing*. Upcoming.}
#'
#' \item{Gotelli, N. J. (2000). Null model analysis of species co-occurrence patterns.
#'  *Ecology, 81(9)*, 2606-2621. <https://doi.org/10.1890/0012-9658(2000)081[2606:NMAOSC]2.0.CO;2>}
#' }
#' @examples
#' \dontrun{
#'
#' my.path <- system.file("extdata", package = "msco")
#' setwd(my.path)
#' my.files <- gtools::mixedsort(list.files(path = my.path, pattern = "*.csv"))
#' my.Archs <- msco::Archetypes(my.files = my.files)
#' my.Archs
#'
#' my.path2 <- system.file("extdata/myCSVs", package = "msco")
#' setwd(my.path2)
#' my.files2 <- gtools::mixedsort(list.files(path = my.path2, pattern = "*.csv"))
#' my.Archs2 <- msco::Archetypes(my.files = my.files2[250:255])
#' my.Archs2
#' }
#' @export
#' @md

Archetypes <- function(my.files, algo = "sim2", metric = "J.occ", nReps = 999){
  mJo.eng(my.files = my.files, algo = algo, metric = metric, nReps = nReps,
            my.Archs = TRUE)$Archs
}

#' A schematic figure of the archetypes
#'
#' A schematic diagram illustrating nine possible archetypes (from the null model test) of the
#'  patterns of species co-occurrences in ecological communities.
#'
#' @return The `Arch_schem` function returns a **schematic diagram** of the archetypes of species
#'  co-occurrence patterns with the following components:
#' \item{**Archetype**}{**Description/Interpretation**}
#' \item{A1}{The joint occupancy value of the observed community matrix (dark solid line) is above the null
#'  model. This means the null hypothesis should be rejected, confirming the presence of a
#'   mechanism of interest being tested (Lagat *et al.,* 202X). It is typical of a community whose
#'   species are positively associated (or aggregated) more often than would be expected by chance.
#'    Such patterns of community structure may arise from a number of ecological mechanisms
#'     including environmental filtering or shared habitat requirements (Cordero and Jackson,
#'      2019).}
#' \item{A2}{A pairwise metric detects a non-random pattern of the community structure,
#'  but when higher order species are considered, a random pattern is produced. This
#'  is typical of a community whose species are aggregated more often than by chance in
#'   sites with few species than in sites with many species (Lagat *et al.,* 202X).}
#' \item{A3}{Species co-occur more often than by chance in sites with few species, but are
#'  segregated more often than by chance in sites with many species, depicting a community
#'   structured by two different community assembly processes (Lagat *et al.,* 202X).}
#' \item{A4}{Using pairwise co-occurrence only fails at detecting patterns of aggregation
#'  for sites with many species, i.e. a type II error (Lagat *et al.,* 202X).}
#' \item{A5}{The test is not statistically significant. This has been ecologically
#'  inferred to mean ecological communities are random and that no community
#'   assembly processes or mechanisms influence their structure (Lagat *et al.,* 202X;
#'    Cordero and Jackson, 2019; Gotelli and Sounding, 2001).}
#' \item{A6}{Using pairwise co-occurrence only fails at detecting patterns of segregation
#'  for sites with many species, i.e. a type II error (Lagat *et al.,* 202X).}
#' \item{A7}{Species are segregated more often than would be expected by chance in sites with
#'  few species, but co-occur more often than would be expected by chance in sites with many
#'   species, depicting a community structured by two different community assembly processes
#'    (Lagat *et al.,* 202X).}
#' \item{A8}{A pairwise metric detects a non-random pattern of the community structure,
#'  but when higher order species are considered, a random pattern is produced. This
#'  is typical of a community whose species are segregated more often than would be expected
#'   by chance in sites with few species than in sites with many species (Lagat *et al.,* 202X).}
#' \item{A9}{The joint occupancy value of the observed community matrix (dark solid line) is below the
#'  null model. This means the null hypothesis should be rejected, confirming the presence of a
#'   mechanism of interest being tested (Lagat *et al.,* 202X). It is typical of a community
#'    structured by inter-specific competition or limiting similarity, though predation might also
#'     generate similar patterns (Hein et al. 2014).}
#' @references
#' \enumerate{
#' \item{Cordero, R.D. and Jackson, D.A. (2019). Species-pair associations, null models, and tests of
#'  mechanisms structuring ecological communities. *Ecosphere* **10.**
#'   <https://doi.org/10.1002/ecs2.2797>}
#'
#'  \item{Gotelli, N. J. and Sounding, E. (2001). Research frontiers in null model analysis. *Glob. Ecol.
#'   Biogeogr.* **10**, 337-343. <https://doi.org/10.1046/j.1466-822X.2001.00249.x>}
#'
#' \item{Hein et al. (2014). Fish introductions reveal the temperature
#'  dependence of species interactions. *Proc. R. Soc. B Biol. Sci.* **281**.
#'   <https://doi.org/10.1098/rspb.2013.2641>}
#' \item{Lagat, V. K., Latombe, G. and Hui, C. (202X). *A multi-species co-occurrence index to
#'  avoid type II errors in null model testing*. Upcoming.}
#' }
#' @note `Arch_schem` is not a generic function which can take in any dataset and give the outputs,
#'  but a path to a schematic diagram saved in this package. A representational figure from empirical,
#'   simulated or any known `.csv` binary data matrices serving as the proxies for species-by-site
#'    presence/absence community data can be accessed with \link[msco]{Jo.plots} function.
#'
#' @examples
#' msco::Arch_schem()
#'
#' @export
#' @md

Arch_schem <- function(){
  schem <- system(paste0('open "', paste0(system.file("ms", package = "msco"), "/Schematic_archetype_figures.pdf"), '"'))
  # magick::image_scale(magick::image_read(system.file(
    # "logos", "arcfig.png", package = "msco")),"597.4x514!")
  return(schem)
}

#' Exponential and power law regression model comparisons
#'
#' Computes the total number of communities with exponential as the best form of joint occupancy
#'  decline than power law and vice versa.
#'
#' @details For a complete list of `AIC` and `Delta_AIC` values of joint occupancy decline regression models
#'  for all communities, see \link[msco]{m.AICs}.
#'
#' The `best.mod2` function is useful when analyzing multiple species-by-site
#'  presence/absence data matrices at once. If one community matrix is analyzed, the `AIC`
#'   output of the function \link[msco]{Jo.eng} should suffice.
#'
#' @param my.files A vector containing names of species-by-site presence/absence `.csv` data matrices.
#'  The data matrices should be saved in the working directory.
#' @return A `table` containig the following components:
#' \item{n}{The number of ecological communities represented by species-by-site
#'  presence/absence `.csv` data matrices.}
#' \item{n.Lwst_AIC}{The number of communities with exponential as the best
#'   form of joint occupancy decline than power law and vice versa.}
#' \item{n.Delta_AIC}{The number of communities whose exponential and power law forms of joint occupancy
#'  decline have `Delta_AIC = 0`, respectively. This number must be equal to `n.Lwst_AIC`.}
#' \item{`%`}{The percentage of `n.Lwst_AIC` (or `n.Delta_AIC`) relative to the total number of
#'  communities (`n`) analyzed.}
#' @references
#' \enumerate{
#' \item{Lagat, V. K., Latombe, G. and Hui, C. (202X). *A multi-species co-occurrence index to
#'  avoid type II errors in null model testing*. Upcoming.}
#'
#' \item{Petrossian, G. A., and Maxfield, M. (2018). An information theory approach to hypothesis testing
#'   in criminological research. *Crime Science*, 7(1), 2.
#'    <https://doi.org/10.1186/s40163-018-0077-5>}
#' }
#' @examples
#' \dontrun{
#'
#' my.path <- system.file("extdata", package = "msco")
#' setwd(my.path)
#' my.files <- gtools::mixedsort(list.files(path = my.path, pattern = "*.csv"))
#' my.best.mod2 <- msco::best.mod2(my.files = my.files)
#' my.best.mod2
#'
#' my.path2 <- system.file("extdata/myCSVs", package = "msco")
#' setwd(my.path2)
#' my.files2 <- gtools::mixedsort(list.files(path = my.path2, pattern = "*.csv"))
#' my2.best.mod2 <- msco::best.mod2(my.files = my.files2[250:255])
#' my2.best.mod2
#' }
#' @export
#' @md

best.mod2 <- function(my.files){
  mJo.eng(my.files, my.best.mod2 = TRUE)$best.mod2
}

#' Exponential, power law and exponential-power law regression model comparisons
#'
#' Computes the total number of communities with exponential or power law or exponential-power law
#'  as the best form of joint occupancy decline among the three (exponential, power law and
#'   exponential-power law) regression models.
#'
#' @details For a complete list of `AIC` and `Delta_AIC` values of joint occupancy decline regression models
#'  for all communities, see \link[msco]{m.AICs}.
#'
#'  The `best.mod3` function is useful when analyzing multiple species-by-site
#'   presence/absence data matrices at once. If one community matrix is analyzed, the `AIC`
#'    output of the function \link[msco]{Jo.eng} should suffice.
#'
#' @param my.files A vector containing names of species-by-site presence/absence `.csv` data matrices.
#'  The data matrices should be saved in the working directory.
#' @return A `table` containig the following components:
#' \item{n}{The number of ecological communities represented by species-by-site
#'  presence/absence `.csv` data matrices.}
#' \item{n.Lwst_AIC}{The number of communities with exponential or power law or
#'   exponential-power law as the best form of joint occupancy decline among the three (exponential,
#'    power law and exponential-power law) regression models.}
#' \item{n.Delta_AIC}{The number of communities whose exponential, power law and exponential-power
#'  law forms of joint occupancy decline, respectively, have `Delta_AIC = 0`. This number must be
#'   equal to `n.Lwst_AIC`.}
#' \item{`%`}{The percentage of `n.Lwst_AIC` (or `n.Delta_AIC`) relative to the total number of
#'  communities (`n`) analyzed.}
#' @references
#' \enumerate{
#' \item{Lagat, V. K., Latombe, G. and Hui, C. (202X). *A multi-species co-occurrence index to
#'  avoid type II errors in null model testing*. Upcoming.}
#'
#' \item{Petrossian, G. A., and Maxfield, M. (2018). An information theory approach to hypothesis testing
#'   in criminological research. *Crime Science*, 7(1), 2.
#'    <https://doi.org/10.1186/s40163-018-0077-5>}
#' }
#' @examples
#' \dontrun{
#'
#' my.path <- system.file("extdata", package = "msco")
#' setwd(my.path)
#' my.files <- gtools::mixedsort(list.files(path = my.path, pattern = "*.csv"))
#' my.best.mod3 <- msco::best.mod3(my.files = my.files)
#' my.best.mod3
#'
#' my.path2 <- system.file("extdata/myCSVs", package = "msco")
#' setwd(my.path2)
#' my.files2 <- gtools::mixedsort(list.files(path = my.path2, pattern = "*.csv"))
#' my2.best.mod3 <- msco::best.mod3(my.files = my.files2[250:255])
#' my2.best.mod3
#' }
#' @export
#' @md

best.mod3 <- function(my.files){
  mJo.eng(my.files, my.best.mod3 = TRUE)$best.mod3
}

#' The Akaike information criterion (AIC) and Delta AIC of joint occupancy decline regression models
#'  for all communities
#'
#' Computes the AIC and Delta AIC of joint occupancy decline regression models for all communities.
#'
#' @details For the number of communities whose joint occupancy decline takes exponential, power
#'  law or exponential-power law forms, see \link[msco]{best.mod2} and \link[msco]{best.mod3}.
#'
#' The `m.AICs` function is useful when analyzing multiple species-by-site
#'  presence/absence data matrices at once. If one community matrix is analyzed, the `AIC`
#'   output of the function \link[msco]{Jo.eng} should suffice.
#'
#' @param my.files A vector containing names of species-by-site presence/absence `.csv` data matrices.
#'  The data matrices should be saved in the working directory.
#' @return A `list` of `data.frame`s containig the following components:
#' \item{df}{The number of parameters in each of the three (exponential, power law and
#'  exponential-power law) joint occupancy decline regression models.}
#' \item{AIC}{The AIC values for each of the three joint occupancy decline regression models.}
#' \item{Delta_AIC3}{The `Delta_AIC` values for each of the three joint occupancy decline regression
#'  models.}
#' \item{Delta_AIC2}{The `Delta_AIC` values for exponential and power law forms of joint occupancy
#'  decline regression models.}
#' @references
#' \enumerate{
#' \item{Lagat, V. K., Latombe, G. and Hui, C. (202X). *A multi-species co-occurrence index to
#'  avoid type II errors in null model testing*. Upcoming.}
#' \item{Petrossian, G. A., and Maxfield, M. (2018). An information theory approach to hypothesis testing
#'   in criminological research. *Crime Science*, 7(1), 2.
#'    <https://doi.org/10.1186/s40163-018-0077-5>}
#'  }
#' @examples
#' \dontrun{
#'
#' my.path <- system.file("extdata", package = "msco")
#' setwd(my.path)
#' my.files <- gtools::mixedsort(list.files(path = my.path, pattern = "*.csv"))
#' my.aic <- msco::m.AICs(my.files = my.files)
#' my.aic
#'
#' my.path2 <- system.file("extdata/myCSVs", package = "msco")
#' setwd(my.path2)
#' my.files2 <- gtools::mixedsort(list.files(path = my.path2, pattern = "*.csv"))
#' my.aic2 <- msco::m.AICs(my.files = my.files2[250:255])
#' my.aic2
#' }
#' @export
#' @md

m.AICs <- function(my.files){
  mJo.eng(my.files, my.AICs = TRUE)$all.AICs
}

#' Parameter estimates of the regression models
#'
#' This function:
#' \itemize{
#' \item{computes the archetypes of the patterns of species co-occurrences in each of the  species-by-site
#'  presence/absence `.csv` data matrices; and}
#' \item{estimates the parameters of:
#'   \enumerate{
#'    \item{**exponential:** **\eqn{J^{\{i\}} = a \times exp(b \times i)}**;}
#'    \item{**power law:** **\eqn{J^{\{i\}} = a \times i^b}**; and}
#'    \item{**exponential-power law:** **\eqn{J^{\{i\}} = a \times exp(b \times i) \times i^c}**}
#'  }}forms of joint occupancy decline, respectively.}
#' @details The `params` function is useful when analyzing multiple species-by-site presence/absence data matrices
#'  at once. If one community matrix is analyzed, the `jo.coeff` and `Archetype` outputs of the function
#'   \link[msco]{Jo.eng} should suffice.
#'
#'
#' @param my.files A vector containing names of species-by-site presence/absence `.csv` data matrices.
#'  The data matrices should be saved in the working directory.
#' @return A `data.frame` consisting of:
#' \item{Arch}{The archetypes of the patterns of species co-occurrences in each of the  species-by-site
#'  presence/absence `.csv` data matrices.}
#' \item{a.ex}{The `a` parameter estimate of the exponential form of joint occupancy decline.}
#' \item{b.ex}{The `b` parameter estimate of the exponential form of joint occupancy decline.}
#' \item{a.pl}{The `a` parameter estimate of the power law form of joint occupancy decline.}
#' \item{b.pl}{The `b` parameter estimate of the power law form of joint occupancy decline.}
#' \item{a.expl}{The `a` parameter estimate of the exponential-power law form of joint occupancy decline.}
#' \item{b.expl}{The `b` parameter estimate of the exponential-power law form of joint occupancy decline.}
#' \item{c.expl}{The `c` parameter estimate of the exponential-power law form of joint occupancy decline.}
#' @references Lagat, V. K., Latombe, G. and Hui, C. (202X). *A multi-species co-occurrence index to
#'  avoid type II errors in null model testing*. Upcoming.
#' @examples
#' \dontrun{
#'
#' my.path <- system.file("extdata", package = "msco")
#' setwd(my.path)
#' my.files <- gtools::mixedsort(list.files(path = my.path, pattern = "*.csv"))
#' my.params <- msco::params(my.files = my.files)
#' my.params
#'
#' my.path2 <- system.file("extdata/myCSVs", package = "msco")
#' setwd(my.path2)
#' my.files2 <- gtools::mixedsort(list.files(path = my.path2, pattern = "*.csv"))
#' my.params2 <- msco::params(my.files = my.files2[250:255])
#' my.params2
#' }
#' @export
#' @md

params <- function(my.files){
  mJo.eng(my.files, my.params = TRUE)$params
}


#' 95% C.I of the parameter estimates of the regression models
#'
#' For every archetype, this function:
#' \itemize{
#' \item{computes the 95% confidence interval of the parameter estimates of:
#'   \enumerate{
#'    \item{**exponential:** **\eqn{J^{\{i\}} = a \times exp(b \times i)}**;}
#'    \item{**power law:** **\eqn{J^{\{i\}} = a \times i^b}**; and}
#'    \item{**exponential-power law:** **\eqn{J^{\{i\}} = a \times exp(b \times i) \times i^c}**}
#'  }}forms of joint occupancy decline, respectively.
#' \item{calculates the exact number of communities under every archetype.}
#' \item{compares the exponential and power law forms of joint occupancy decline using their
#'  AIC values.}
#' \item{contrasts the three parametric models (exponential, power law and exponential-power law)
#'    of joint occupancy decline using their AIC values.}
#'  }
#' @details The `params_c.i` function is useful when analyzing multiple species-by-site presence/absence
#'  data matrices at once. If one community matrix is analyzed, the `jo.coeff` and `Archetype` outputs
#'   of the function \link[msco]{Jo.eng} should suffice.
#'
#' @param my.files A vector containing names of species-by-site presence/absence `.csv` data matrices.
#'  The data matrices should be saved in the working directory.
#' @return A `data.frame` consisting of:
#' \item{Arch}{The archetypes of the patterns of species co-occurrences in each of the  species-by-site
#'  presence/absence `.csv` data matrices.}
#' \item{n}{The number of communities under every archetype.}
#' \item{`ex_%`}{The percentages of the number of communities (under every archetype) where
#'   exponential form of joint occupancy decline fitted better than power law.}
#' \item{a.ex}{The 95% closed confidence interval of the `a` parameter estimates of the exponential
#'  form of joint occupancy decline, under every archetype.}
#' \item{b.ex}{The 95% closed confidence interval of the `b` parameter estimates of the exponential
#'  form of joint occupancy decline, under every archetype.}
#' \item{`p.l_%`}{The percentages of the number of communities (under every archetype) where
#'   power law form of joint occupancy decline fitted better than exponential.}
#' \item{a.pl}{The 95% closed confidence interval of the `a` parameter estimates of the power law
#'  form of joint occupancy decline, under every archetype.}
#' \item{b.pl}{The 95% closed confidence interval of the `b` parameter estimates of the power law
#'  form of joint occupancy decline, under every archetype.}
#' \item{`ex.pl_%`}{The percentages of the number of communities (under every archetype) where exponential-power
#'    law form of joint occupancy decline fitted better than both the exponential and power law forms.}
#' \item{a.expl}{The 95% closed confidence interval of the `a` parameter estimates of the exponential-power law
#'  form of joint occupancy decline, under every archetype.}
#' \item{b.expl}{The 95% closed confidence interval of the `b` parameter estimates of the exponential-power law
#'  form of joint occupancy decline, under every archetype.}
#' \item{c.expl}{The 95% closed confidence interval of the `c` parameter estimates of the exponential-power law
#'  form of joint occupancy decline, under every archetype.}
#' @references Lagat, V. K., Latombe, G. and Hui, C. (202X). *A multi-species co-occurrence index to
#'  avoid type II errors in null model testing*. Upcoming.
#' @examples
#' \dontrun{
#'
#' my.path <- system.file("extdata", package = "msco")
#' setwd(my.path)
#' my.files <- gtools::mixedsort(list.files(path = my.path, pattern = "*.csv"))
#' my.params_c.i <- msco::params_c.i(my.files = my.files)
#' my.params_c.i
#'
#' my.path2 <- system.file("extdata/myCSVs", package = "msco")
#' setwd(my.path2)
#' my.files2 <- gtools::mixedsort(list.files(path = my.path2, pattern = "*.csv"))
#' my.params_c.i2 <- msco::params_c.i(my.files = my.files2[250:255])
#' my.params_c.i2
#' }
#' @export
#' @md

params_c.i <- function(my.files){
  mJo.eng(my.files, my.params_c.i = TRUE)$params_c.i
}

#' Joint occupancy parametric and null model plots
#'
#' Plots the null model and joint occupancy decline with order (number of species)
#' and fits the decline to exponential, power law and exponential-power law
#' parametric models, respectively.
#'
#' @param jo_Obj A joint occupancy model object of `class:` `jo_Obj`
#' @return Produces a figure consisting of the following plots:
#'  \item{(a)}{Joint occupancy decline.}
#'  \item{(b)}{Exponential regression of the joint occupancy decline.}
#'  \item{(c)}{Power law regression of the joint occupancy decline.}
#'  \item{(d)}{Exponential-power law regression of the joint occupancy decline.}
#'  \item{(e)}{Null model test.}
#' @details This function provides a clear visualization of the forms of joint
#'   occupancy decline and null model test. It offers information on:
#'  \itemize{
#'  \item{the outcomes of the null model test (through the appended
#'    archetype value on (e) plot) and}
#'  \item{the comparisons between the forms of joint occupancy decline (through
#'      the affixed `AIC` and `rsq` values on (b), (c) and (d) plots,
#'       respectively)}.
#'   }
#' @references Lagat, V. K., Latombe, G. and Hui, C. (202X). *A multi-species co-occurrence
#'  index to avoid type II errors in null model testing*. Upcoming.
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

#' Joint occupancy parametric and null model plots for multiple communities
#'
#' For multiple communities, this function plots the null model and joint
#'  occupancy decline with order (number of species) and fits the decline
#'   to exponential, power law and exponential-power law parametric models,
#'    respectively.
#'
#' @param my.files A vector containing names of species-by-site presence/absence `.csv` data matrices.
#'  The data matrices should be saved in the working directory.
#' @param algo Simulation algorithm used. The possible options to choose from are: `sim1`,
#'  `sim2`, `sim3`, `sim4`, `sim5`, `sim6`, `sim7`, `sim8`, and `sim9`, all from
#'   Gotelli (2000). `sim2` is highly recommended (see Lagat *et al.,* 202X).
#' @param metric Metric used to quantify the patterns in `s.data`. It has to be multi-species
#'  co-occurrence index (see \link[msco]{j.occ}).
#' @param nReps Number of simulations used in the null model test.
#' @return Produces a `.pdf` file with multiple figures each consisting of the
#'  following plots:
#'  \item{(a)}{as for \link[msco]{Jo.plots}}
#'  \item{(b)}{as for \link[msco]{Jo.plots}}
#'  \item{(c)}{as for \link[msco]{Jo.plots}}
#'  \item{(d)}{as for \link[msco]{Jo.plots}}
#'  \item{(e)}{as for \link[msco]{Jo.plots}}
#'
#' @details This function provides a clear visualization of the forms of joint
#'   occupancy decline and null model test, for multiple communities. For more
#'    details on individual community joint occupancy parametric and null model
#'     plots, see \link[msco]{Jo.plots}.
#'
#'  The `m.Jo.plots` function is useful when analyzing multiple species-by-site
#'   presence/absence data matrices at once. If one community matrix is analyzed,
#'   the output of the function \link[msco]{Jo.plots} should suffice.
#'
#' @references
#' \enumerate{
#' \item{Lagat, V. K., Latombe, G. and Hui, C. (202X). *A multi-species co-occurrence
#'  index to avoid type II errors in null model testing*. Upcoming.}
#'
#' \item{Gotelli, N. J. (2000). Null model analysis of species co-occurrence patterns.
#'  *Ecology, 81(9)*, 2606-2621. <https://doi.org/10.1890/0012-9658(2000)081[2606:NMAOSC]2.0.CO;2>}
#'  }
#'
#' @examples
#' \dontrun{
#'
#' my.path <- system.file("extdata", package = "msco")
#' setwd(my.path)
#' my.files <- gtools::mixedsort(list.files(path = my.path, pattern = "*.csv"))
#' m.plots <- msco::m.Jo.plots(my.files = my.files)
#' m.plots
#'
#' ## Close the open pdf reader
#' my.path2 <- system.file("extdata/myCSVs", package = "msco")
#' setwd(my.path2)
#' my.files2 <- gtools::mixedsort(list.files(path = my.path2, pattern = "*.csv"))
#' m.plots2 <- msco::m.Jo.plots(my.files = my.files2[250:255])
#' m.plots2
#' }
#' @export
#' @md

m.Jo.plots <- function(my.files, algo = "sim2", metric = "J.occ", nReps = 999){
  mJo.eng(my.files,  algo = algo, metric = metric, nReps = nReps, m.Jo.plots = TRUE)$m.Jo.plots
}

#' Robustness of joint occupancy decline regression models
#'
#' Determines the robustness of the exponential, power law and exponential-power law
#'  forms of joint occupancy decline (for multiple communities) by computing the square
#'   of Pearson's correlation coefficient (\eqn{r}) between the joint occupancy values
#'    of the observed data and predicted data, for all orders of species.
#'
#' @details The `rsq` function is useful when analyzing multiple species-by-site
#'  presence/absence data matrices at once. If one community matrix is analyzed, the `r2`
#'   output of the function \link[msco]{Jo.eng} should suffice.
#'
#' @param my.files A vector containing names of species-by-site presence/absence `.csv` data matrices.
#'  The data matrices should be saved in the working directory.
#' @return A `list` of `data.frame`s containig the following components:
#' \item{`rsq.ex`}{\eqn{r^2} for the exponential form of joint occupancy decline.}
#' \item{rsq.pl}{\eqn{r^2} for the power law form of joint occupancy decline.}
#' \item{rsq.ex.pl}{\eqn{r^2} for the exponential-power law form of joint occupancy decline.}
#' @references
#' \enumerate{
#'
#' \item{Lagat, V. K., Latombe, G. and Hui, C. (202X). *A multi-species co-occurrence
#'  index to avoid type II errors in null model testing*. Upcoming.}
#'
#' \item{Pearson, K. (1895) VII. Note on regression and inheritance in the
#'  case of two parents. *proceedings of the royal society of London,* **58**:240-242.
#'   <https://doi.org/10.1098/rspl.1895.0041>}
#' }
#' @seealso \link[msco]{rsq.s}
#' @examples
#' \dontrun{
#'
#' my.path <- system.file("extdata", package = "msco")
#' setwd(my.path)
#' my.files <- gtools::mixedsort(list.files(path = my.path, pattern = "*.csv"))
#' my.r2 <- msco::rsq(my.files = my.files)
#' my.r2
#'
#' my.path2 <- system.file("extdata/myCSVs", package = "msco")
#' setwd(my.path2)
#' my.files2 <- gtools::mixedsort(list.files(path = my.path2, pattern = "*.csv"))
#' my2.r2 <- msco::rsq(my.files = my.files2[250:255])
#' my2.r2
#' }
#' @export
#' @md

rsq <- function(my.files){
  mJo.eng(my.files, my.r2 = TRUE)$r2
}

#' Summary: Robustness of joint occupancy decline regression models
#'
#' Gives a summary of the total number of communities (under each and for all
#'  archetypes) whose forms of joint occupancy decline have \eqn{r^2 > 0.95}.
#' @details The `rsq.s` function is useful when analyzing multiple species-by-site
#'  presence/absence data matrices at once. If one community matrix is analyzed, the `r2`
#'   output of the function \link[msco]{Jo.eng} should suffice.
#'
#' @param my.files A vector containing names of species-by-site presence/absence `.csv` data matrices.
#'  The data matrices should be saved in the working directory.
#' @return A `list` containig the following components:
#'
#'  $`rsq.per.Archs`
#'
#' \item{`Archs`}{Archetypes of the patterns of species co-occurrences in each
#'  of the species-by-site presence/absence .csv data matrices.}
#' \item{`n.a`}{Number of communities under each archetype.}
#' \item{`rsq.ex`}{Number of communities under each archetype whose exponential
#'  forms of joint occupancy decline have \eqn{r^2 > 0.95}.}
#' \item{`rsq.pl`}{Number of communities under each archetype whose power
#'  law forms of joint occupancy decline have \eqn{r^2 > 0.95}.}
#' \item{`rsq.ex-pl`}{Number of communities under each archetype whose
#'  exponential-power law forms of joint occupancy decline have \eqn{r^2 > 0.95}.}
#'
#'  $`rsq.all.Communities`
#'
#' \item{n}{Number of all communities analyzed}
#' \item{`ex`}{Number of communities whose exponential forms of joint occupancy
#'  decline have \eqn{r^2 > 0.95}}
#' \item{`pl`}{Number of communities whose power law forms of joint occupancy
#'  decline have \eqn{r^2 > 0.95}}
#' \item{`ex.pl`}{Number of communities whose exponential-power law forms
#'  of joint occupancy decline have \eqn{r^2 > 0.95}}
#'
#' @references
#' \enumerate{
#'
#' \item{Lagat, V. K., Latombe, G. and Hui, C. (202X). *A multi-species co-occurrence
#'  index to avoid type II errors in null model testing*. Upcoming.}
#'
#' \item{Pearson, K. (1895) VII. Note on regression and inheritance in the
#'  case of two parents. *proceedings of the royal society of London,* **58**:240-242.
#'   <https://doi.org/10.1098/rspl.1895.0041>}
#' }
#' @seealso \link[msco]{rsq}
#' @examples
#' \dontrun{
#'
#'my.path <- system.file("extdata", package = "msco")
#' setwd(my.path)
#' my.files <- gtools::mixedsort(list.files(path = my.path, pattern = "*.csv"))
#' my.rsq.s <- msco::rsq.s(my.files = my.files)
#' my.rsq.s
#'
#' my.path2 <- system.file("extdata/myCSVs", package = "msco")
#' setwd(my.path2)
#' my.files2 <- gtools::mixedsort(list.files(path = my.path2, pattern = "*.csv"))
#' my2.rsq.s <- msco::rsq.s(my.files = my.files2[250:255])
#' my2.rsq.s
#' }
#' @export
#' @md

rsq.s <- function(my.files){
  mJo.eng(my.files, my.r2.s = TRUE)$r2.s
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
#' @return Returns a `list` with the following outputs:
#' \item{jo.val}{Joint occupancy value.}
#' \item{jo.sd}{The standard deviation of `jo.val`.}
#' @references Lagat, V. K., Latombe, G. and Hui, C. (202X). *A multi-species co-occurrence
#'  index to avoid type II errors in null model testing*. Upcoming.
#' @examples
#' ex.data <- read.csv(system.file("extdata", "274.csv", package = "msco"))
#' jo <- msco::j.occ(ex.data, order = 3)
#' jo
#'
#' ex.data2 <- read.csv(system.file("extdata", "65.csv", package = "msco"))
#' jo2 <- msco::j.occ(ex.data2, order = 3)
#' jo2
#' @export
#' @md

j.occ<-function(s.data, order){

  s.data <- as.matrix(s.data)
  richness <- colSums(s.data)
  p <- exp(lchoose(richness, order) - lchoose(nrow(s.data),order))
  jo.val <- sum(p)

  similarity_mat <- t(s.data) %*% s.data
  covmat <- exp(lchoose(similarity_mat, order) - lchoose(nrow(s.data),
                                                         order))
  for (j in 1:ncol(s.data)) {
    for (k in 1:ncol(s.data)) {
      covmat[j, k] <- covmat[j, k] - p[j] * p[k]
    }
  }
  jo.var <- choose(nrow(s.data), order)/
    (choose(nrow(s.data), order) - 1) * sum(covmat)

  jo.sd <- sqrt(jo.var)

  jo.index <- list()
  jo.index$jo.val <- jo.val
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
#' @return Returns a `list` with the following outputs:
#' \item{jo.vals}{A vector of joint occupancy values for a range number of species (in `orders`).}
#' \item{jo.sds}{A vector of standard deviations of `jo.vals`.}
#' @references Lagat, V. K., Latombe, G. and Hui, C. (202X). *A multi-species co-occurrence
#'  index to avoid type II errors in null model testing*. Upcoming.
#' @examples
#' ex.data <- read.csv(system.file("extdata", "274.csv", package = "msco"))
#' jos <- msco::j.occs(ex.data, orders = 1:nrow(ex.data))
#' jos
#'
#' ex.data2 <- read.csv(system.file("extdata", "65.csv", package = "msco"))
#' jos2 <- msco::j.occs(ex.data2, orders = 1:nrow(ex.data))
#' jos2
#' @export
#' @md

j.occs<-function(s.data, orders = 1:nrow(s.data)){

  jo=0
  SDs=0
  for (i in orders) {
    jo[i]=j.occ(s.data, order = i)$jo.val
    SDs[i]=j.occ(s.data, order = i)$jo.sd
  }
  jo.inds <- list()
  jo.inds$jo.vals <- jo
  jo.inds$jo.sds <- SDs
  return(jo.inds)
}

#' Results on joint occupancy index (presented in Lagat et al., 202Xa)
#'
#' This function is not for general use but avails the results presented
#'  in Lagat *et al*. (202Xa).
#'
#' @return Returns all the results presented in Lagat *et al*. (202Xa). To replicate these results,
#'  execute the following code:
#'
#'  ```
#'    RNGkind(sample.kind = "Rejection")
#'    set.seed(14)
#'    my.path <- system.file("extdata/myCSVs", package = "msco")
#'    setwd(my.path)
#'    my.files <- gtools::mixedsort(list.files(path = my.path, pattern = "*.csv"))
#'    Lag.res <- msco::mJo.eng(my.files,
#'                     m.Jo.plots = TRUE,
#'                     my.Archs = FALSE,
#'                     my.AICs = FALSE,
#'                     my.params = FALSE,
#'                     my.r2 = FALSE,
#'                     my.r2.s = TRUE,
#'                     my.best.mod2 = TRUE,
#'                     my.best.mod3 = TRUE,
#'                     my.params_c.i = TRUE
#'                    );Lag.res
#'
#'  ```
#'
#'  **Caveat:** The above code can take approximately 10 minutes to execute. It took 10.39014 minutes to run
#'   (and output results) on a 64 bit system with 8 GB RAM and 3.60 GHz CPU.
#'
#'  * __Fig. 3__ can be replicated using:
#'  ```
#'   RNGkind(sample.kind = "Rejection")
#'   set.seed(14)
#'   path <- system.file("ms", package = "msco")
#'   grDevices::pdf(file = paste0(path, "/real.arch.plots.pdf"),
#'    paper="a4r", height = 8.27, width = 11.69)
#'   msco:::nullmod_archs()
#'   grDevices::dev.off()
#'   system(paste0('open "', paste0(path, "/real.arch.plots.pdf"), '"'))
#'  ```
#'
#' @references Lagat, V. K., Latombe, G. and Hui, C. (202Xa). *A multi-species co-occurrence
#'  index to avoid type II errors in null model testing*. Upcoming.
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
  res <- readRDS(system.file("ms", "m.Jo.res.RDS", package = "msco"))
  # saveRDS(j.en, file = paste0(system.file("ms", package="msco"), "/jo.res.RDS"))
  # res.plo <- Biobase::openPDF(system.file("ms", "jo.plots.pdf", package = "msco"))
  system(paste0('open "', paste0(system.file("ms", package = "msco"), "/real.arch.plots.pdf"), '"'))
  system(paste0('open "', paste0(system.file("ms", package = "msco"), "/m.Jo.plots.pdf"), '"'))

  return(res)
}

#' Results on generalised B-spline modelling (presented in Lagat et al., 202Xb)
#'
#' This function is not for general use but avails the generalised B-spline
#'  model results presented in Lagat *et al*. (202Xb).
#'
#' @return Returns all the results presented in Lagat *et al*. (202Xb). To replicate
#'
#' * __Figs. 1__, __3__, __4__, __5__, and __Tables 1__ and __S1__, execute the following code:
#'
#'
#'   ```
#'
#'    my.path <- system.file("extdata/gsmdat", package = "msco")
#'    setwd(my.path)
#'    s.data <- get(load("s.data.csv")) # Species-by-site matrix
#'    t.data <- get(load("t.data.csv")) # Species-by-trait matrix
#'    p.d.mat <- get(load("p.d.mat.csv")) # Species-by-species phylogenetic distance matrix
#'    RNGkind(sample.kind = "Rejection")
#'    set.seed(0)
#'    gb.res <- msco::gbsm_m.orders(s.data,
#'                t.data,
#'                p.d.mat,
#'                metric="Simpson_eqn",
#'                orders=c(2:5, 8, 10, 15),
#'                d.f=4,
#'                degree=3,
#'                n=1000,
#'                k=5,
#'                p=0.8,
#'                type="k-fold",
#'                scat.plots=TRUE,
#'                response.curves=TRUE,
#'                j.occs.distrbn=TRUE,
#'                mp.plots=TRUE,
#'                start=seq(-0.1, 0, length.out=(ncol(t.data)+2)*4+1)
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
#'        s.data <- get(load("s.data.csv")) ## Species-by-site matrix
#'        taxa <- get(load("taxa.levels.csv")) ## Species taxa
#'        my.phylo.plot <- msco::s.phylo(s.data,
#'                                  database = "ncbi",
#'                                  obs.taxa=FALSE,
#'                                  taxa.levels = taxa,
#'                                  Obs.data=FALSE,
#'                                  phy.d.mat=FALSE,
#'                                  phylo.plot = TRUE)
#'
#'      ```
#'
#'      + __Fig. 2__:
#'      ```
#'        my.path <- system.file("extdata/gsmdat", package = "msco")
#'        setwd(my.path)
#'        s.data <- get(load("s.data.csv")) ## Species-by-site matrix
#'        t.data <- get(load("t.data.csv")) ## Species-by-trait matrix
#'        p.d.mat <- get(load("p.d.mat.csv")) ## Species-by-species pgylogenetic distance matrix
#'        my.gbsm <- msco::gbsm(s.data,
#'                          t.data,
#'                          p.d.mat,
#'                          metric = "Simpson_eqn",
#'                          d.f=4,
#'                          order.jo=3,
#'                          degree=3,
#'                          n=1000,
#'                          b.plots=TRUE,
#'                          bsplines="single",
#'                          scat.plot=FALSE,
#'                          response.curves=FALSE,
#'                          leg=1,
#'                          start=seq(-0.1, 0, length.out=(ncol(t.data)+2)*4+1)
#'                        )
#'
#'      ```
#'      + __Fig. S2__:
#'      ```
#'        my.path <- system.file("extdata/gsmdat", package = "msco")
#'        setwd(my.path)
#'        s.data <- get(load("s.data.csv")) ## Species-by-site matrix
#'        t.data <- get(load("t.data.csv")) ## Species-by-trait matrix
#'        p.d.mat <- get(load("p.d.mat.csv")) ## Species-by-species phylogenetic distance matrix
#'        RNGkind(sample.kind = "Rejection")
#'        set.seed(0)
#'        pe <- msco::pred.error.bands(s.data,
#'                            t.data,
#'                            p.d.mat,
#'                            metric="Simpson_eqn",
#'                            d.f=4,
#'                            simm=10,
#'                            orders = c(2:5, 8, 10, 15),
#'                            degree=3,
#'                            n=1000,
#'                            start=seq(-0.1, 0, length.out=(ncol(t.data)+2)*4+1)
#'                          )
#'      ```
#'
#'  **Caveat:** The above codes can collectively take approximately 7 minutes to execute
#'   (with prediction uncertainty plot taking 6 minutes alone). It took 7.3895 minutes
#'    to run (and output results) on a 64 bit system with 8 GB RAM and 3.60 GHz CPU.
#'
#' @references Lagat, V. K., Latombe, G. and Hui, C. (202Xb). *Dissecting the effects of
#'   neutral encounter versus functional traits on multi-order species interactions
#'    and co-occurrence with generalised B-spline modelling*. Upcoming.
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
  # Biobase::openPDF(system.file("ms", "gbsm.plots.pdf", package = "msco"))
  system(paste0('open "', paste0(system.file("ms", package = "msco"), "/gbsm.plots.pdf"), '"'))
  # Biobase::openPDF(system.file("ms", "Cluster.Dendrogram.pdf", package = "msco"))
  # Biobase::openPDF(system.file("ms", "B-splines.curves.pdf", package = "msco"))
  # Biobase::openPDF(system.file("ms", "pred.error.bands.pdf", package = "msco"))
  # Biobase::openPDF(system.file("ms", "gbsm.plots.pdf", package = "msco"))
  # saveRDS(gb.res, file = paste0(system.file("ms", package="msco"), "/gbsm.res.RDS"))

  return(gres)
}

#' Results on `msco` illustration (presented in Lagat et al., 202Xc)
#'
#' This function is not for general use but avails the `msco` R package
#'  illustration presented in Lagat *et al*. (202Xc).
#'
#' @return Returns all the results presented in Lagat *et al*. (202Xc). To replicate
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
#'                algo="sim2",
#'                metric = "j.occ",
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
#'                lab=FALSE,
#'                leg=FALSE,
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
#'    j.en$nmod_stats ## Table 2
#'    grDevices::pdf(file = paste0(system.file("ms", package = "msco"),
#'     "/aJo.plots.pdf"), paper="a4r", height = 8.27, width = 11.69)
#'    j.en$all.plots
#'    grDevices::dev.off()
#'    system(paste0('open "', paste0(system.file("ms", package = "msco"), ## Fig. 2
#'     "/aJo.plots.pdf"), '"'))
#'
#'   ```
#'
#' * __Fig. 3__, execute the following code:
#'
#'   ```
#'    RNGkind(sample.kind = "Rejection")
#'    set.seed(14)
#'    grDevices::pdf(file = paste0(system.file("ms", package = "msco"),
#'     "/real.arch.plots2.pdf"), paper="a4r", height = 8.27, width = 11.69)
#'    msco:::nullmod_archs2()
#'    grDevices::dev.off()
#'    system(paste0('open "', paste0(system.file("ms", package = "msco"),
#'     "/real.arch.plots2.pdf"), '"'))
#'
#'   ```
#'
#' * __Fig. 5__, execute the following code:
#'
#'   ```
#'
#'    my.path <- system.file("extdata/gsmdat", package = "msco")
#'    setwd(my.path)
#'    s.data <- get(load("s.data.csv")) # Species-by-site matrix
#'    t.data <- get(load("t.data.csv")) # Species-by-trait matrix
#'    p.d.mat <- get(load("p.d.mat.csv")) # Species-by-species phylogenetic distance matrix
#'    RNGkind(sample.kind = "Rejection")
#'    set.seed(0)
#'    gb.res <- msco::gbsm_m.orders(s.data,
#'                t.data,
#'                p.d.mat,
#'                metric="Simpson_eqn",
#'                orders=c(3:5, 8, 10, 15, 20),
#'                d.f=4,
#'                degree=3,
#'                n=1000,
#'                k=5,
#'                p=0.8,
#'                type="k-fold",
#'                scat.plots=FALSE,
#'                response.curves=TRUE,
#'                j.occs.distrbn=FALSE,
#'                mp.plots=FALSE,
#'                start=seq(-0.1, 0, length.out=(ncol(t.data)+2)*4+1)
#'              )
#'
#'
#'   ```
#' @references
#' \enumerate{
#' \item{Lagat, V. K., Latombe, G. and Hui, C. (202Xa). *A multi-species co-occurrence
#'  index to avoid type II errors in null model testing*. Upcoming.}
#'
#'  \item{Lagat, V. K., Latombe, G. and Hui, C. (202Xb). *Dissecting the effects of
#'   neutral encounter versus functional traits on multi-order species interactions
#'    and co-occurrence with generalised B-spline modelling*. Upcoming.}
#'
#'   \item{Lagat, V. K., Latombe, G. and Hui, C. (202Xc). `msco:` *an R software
#'    package for null model testing of multi-species interactions and  interference,
#'     and analysis of their drivers*. Upcoming.}
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
  system(paste0('open "', paste0(system.file("ms", package = "msco"), "/msco.illus.plots.pdf"), '"'))
  # saveRDS(j.en, file = paste0(system.file("ms", package="msco"), "/msco.res.RDS"))

  return(mres)
}


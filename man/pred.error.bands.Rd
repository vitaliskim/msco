% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pred.error.bands.R
\name{pred.error.bands}
\alias{pred.error.bands}
\title{Prediction uncertainty}
\usage{
pred.error.bands(
  s.data,
  t.data,
  p.d.mat,
  metric = "Simpson_eqn",
  gbsm.model,
  d.f = 4,
  simm = 10,
  orders,
  degree = 3,
  n = 1000,
  max.vif = 40,
  max.vif2 = 30,
  start.range = c(-0.1, 0)
)
}
\arguments{
\item{s.data}{A species-by-site presence/absence \code{data.frame} with entries indicating
occurrence (1) and non-occurrence (0) of species in a site.}

\item{t.data}{A \code{data.frame} with traits as columns and species as rows. The species must be the same as
in \code{s.data}.}

\item{p.d.mat}{A symmetric \code{matrix} with \code{dimnames} as species and entries indicating the
phylogenetic distance between any two of them (species).}

\item{metric}{As for \link[msco]{gbsm_m.orders}.}

\item{gbsm.model}{As for \link[msco]{gbsm_m.orders}.}

\item{d.f}{As for \link[msco]{gbsm_m.orders}.}

\item{simm}{Number of Monte Carlo simulations performed}

\item{orders}{As for \link[msco]{gbsm_m.orders}}

\item{degree}{As for \link[msco]{gbsm_m.orders}.}

\item{n}{As for \link[msco]{gbsm_m.orders}.}

\item{max.vif}{As for \link[msco]{gbsm}.}

\item{max.vif2}{As for \link[msco]{gbsm}.}

\item{start.range}{As for \link[msco]{gbsm_m.orders}.}
}
\value{
\code{pred.error.bands} function returns:
\item{\code{predictors}}{a \code{data.frame} of predictors}
\item{\code{responses}}{a \code{data.frame} of response values of predictors}
\item{\code{responses.sim_stats}}{a \code{data.frame} of the reponses' mean and standard deviation
(from \code{simm} replicates), and}
\itemize{
\item the response curves with prediction error bands for all orders of joint occupancy
}
}
\description{
This function plots the response curves showing the effect of the predictors (i.e. trait-based and
neutral forces) on joint occupancy as the response variable, with prediction error bands (as the
standard deviation from the mean of the response variable) for all orders of joint occupancy.
}
\examples{
\dontrun{
 my.path <- system.file("extdata/gsmdat", package = "msco")
 setwd(my.path)
 s.data <- get(load("s.data.csv")) ## Species-by-site matrix
 t.data <- get(load("t.data.csv")) ## Species-by-Trait matrix
 p.d.mat <- get(load("p.d.mat.csv")) ## Species-by-species phylogenetic distance matrix

 RNGkind(sample.kind = "Rejection")
 set.seed(1)
 pe <- msco::pred.error.bands(s.data, t.data, p.d.mat, metric="Simpson_eqn", d.f=4, simm=10,
  orders = c(2:5, 8, 10, 15), degree=3, n=1000, gbsm.model, start.range=c(-0.2, 0))

 pe$predictors$`order 2`
 pe$responses$`order 2`
 pe$responses.sim_stats$`order 2`

 pe$predictors$`order 3`
 pe$responses$`order 3`
 pe$responses.sim_stats$`order 3`

 pe$predictors$`order 10`
 pe$responses$`order 10`
 pe$responses.sim_stats$`order 10`

 }

}
\references{
\enumerate{
\item{Lagat, V. K., Latombe, G. and Hui, C. (2021a). \emph{A multi-species co-occurrence
index to avoid type II errors in null model testing}. DOI: \verb{<To be added>}.}

\item{Lagat, V. K., Latombe, G. and Hui, C. (2021b). \emph{Dissecting the effects of random
encounter versus functional trait mismatching on multi-species co-occurrence and
interference with generalised B-spline modelling}. DOI: \verb{<To be added>}.}

}
}

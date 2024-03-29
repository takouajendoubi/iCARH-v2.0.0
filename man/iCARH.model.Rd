% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iCARH.model.R
\name{iCARH.model}
\alias{iCARH.model}
\title{Runs the integrative CAR Horseshoe model}
\usage{
iCARH.model(X, Y, drug, groups = NULL, pathways, tau = 1.2,
  NA_value = -99999, init = T, ...)
}
\arguments{
\item{X}{the metabolomics time-course data}

\item{Y}{the additional omic time-course data}

\item{drug}{treatment effect (NA values not allowed in drug).}

\item{groups}{grouping vector (discrete)}

\item{pathways}{pathway adjacency matrices}

\item{tau}{global sparsity parameter \eqn{\tau} as in Jendoubi, T., & Ebbels, T. (2018)}

\item{NA_value}{NA values are incompatible with stan.
NAs will be replaced by NA_value and will be inferred (only for X and Y data).}

\item{init}{If \code{TRUE} use iCARH provided initialization function. Passed to Stan otherwise. Please see Stan manual 
on \code{init} possible values.}

\item{...}{additional stan parameters}
}
\value{
stan object
}
\description{
Infers treatment effects, association with heterogeneous omic variables, pathway perturbation
among other parameters (e.g. time dependence). Regression coefficients (beta parameter) are initialized
using a univariate regression ignoring time and metabolite dependence.
}
\examples{
data.sim = iCARH.simulate(4, 8, 10, 2, 2, path.probs=0.3, Zgroupeff=c(0,4),
beta.val=c(1,-1,0.5, -0.5))
XX = data.sim$XX
Y = data.sim$Y
Z = data.sim$Z
pathways = data.sim$pathways
\donttest{
rstan_options(auto_write = TRUE)
options(mc.cores = 2)
fit = iCARH.model(XX, Y, Z, pathways, control = list(adapt_delta = 0.99, max_treedepth=10),
iter = 2, chains = 2)}

}

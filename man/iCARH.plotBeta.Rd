% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/postprocessing.R
\name{iCARH.plotBeta}
\alias{iCARH.plotBeta}
\alias{iCARH.plotTreatmentEffect}
\alias{iCARH.plotPathwayPerturbation}
\alias{iCARH.plotDataImputation}
\alias{iCARH.checkRhats}
\alias{iCARH.checkNormality}
\alias{iCARH.waic}
\alias{iCARH.mad}
\title{Postprocess and plot model parameters}
\usage{
iCARH.plotBeta(fit, indx = TRUE, indy = TRUE)

iCARH.plotTreatmentEffect(fit, indx = TRUE)

iCARH.plotPathwayPerturbation(fit, path.names, indpath = TRUE)

iCARH.plotDataImputation(fit, indx = T, indy = T, plotx = T,
  ploty = T, ...)

iCARH.checkRhats(fit)

iCARH.checkNormality(fit)

iCARH.waic(fit)

iCARH.mad(fit)
}
\arguments{
\item{fit}{object returned by iCARH.model}

\item{indx}{vector to specify X variables to plot. Selects all variables of X by default.}

\item{indy}{vector to specify Y variables to plot. Selects all variables of Y by default.}

\item{path.names}{pathway names}

\item{indpath}{vector to specify pathways to plot. Selects all pathways by default.}

\item{plotx}{plot X data imputation?}

\item{ploty}{plot Y data imputation?}

\item{...}{passed to ggplot2::geom_violin}
}
\value{
the \code{iCARH.plot[*]} functions return a ggplot graph object. \code{iCARH.checkNormality} returns the normalized data.
\code{iCARH.waic} and \code{iCARH.mad} return corresponding waic (scalar) and mad (vector of \eqn{J*(J+1)/2}) values. 
\code{iCARH.checkRhats} checks model convergence.
}
\description{
Group of functions to postprocess and plot model parameters of interest, compute WAIC
(Watanabe-Akaike Information Criterion) and MADs (Mean Absolute Deviation) for posterior predictive checks
and check normality assumptions.
}
\section{Functions}{
\itemize{
\item \code{iCARH.plotBeta}: Plot boxplots of posterior densities of \eqn{\beta} coefficients.

\item \code{iCARH.plotTreatmentEffect}: Plot boxplots of posterior densities of treatment effect coefficients.

\item \code{iCARH.plotPathwayPerturbation}: Plot posterior densities of pathway perturbation parameters

\item \code{iCARH.plotDataImputation}: Plot imputed values

\item \code{iCARH.checkRhats}: check model convergence and return Rhat coefficients

\item \code{iCARH.checkNormality}: Check normality assumptions. Returns normalized data
and performs quantile-quantile plot

\item \code{iCARH.waic}: Compute Watanabe-Akaike Information Criterion (WAIC)

\item \code{iCARH.mad}: Compute MADs (Mean Absolute Deviation) between true covariance matrix
and inferred covariance matrix for posterior predictive checks
}}

\examples{
data.sim = iCARH.simulate(4, 10, 14, 8, 2, path.probs=0.3, Zgroupeff=c(0,4),
beta.val=c(1,-1,0.5, -0.5))
XX = data.sim$XX
Y = data.sim$Y
Z = data.sim$Z
pathways = data.sim$pathways
\donttest{
rstan_options(auto_write = TRUE)
options(mc.cores = 2)
fit = iCARH.model(XX, Y, Z, pathways, control = list(adapt_delta = 0.99, max_treedepth=10),
iter = 2, chains = 2)
if(!is.null(fit))
gplot = iCARH.plotBeta(fit, indx=1:3, indy=1:2)}

}

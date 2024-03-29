% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{noisy_data}
\alias{noisy_data}
\title{Temperature exceedance probabilities and observations}
\format{
A dataframe containing 8 variables:
\describe{
  \item{Obs}{Values equal to 0 or 1 depending on whether the observed temperature
  does (\code{Obs = 1}) or does not (\code{Obs = 0}) exceed the threshold \code{th}.}
  \item{Ens}{Probability forecast issued by the COSMO-E ensemble prediction system
  when predicting whether the observed temperature exceeds \code{th} (i.e. \code{Obs = 1}).}
  \item{Clim}{Probability forecast issued by a climatological forecast distribution
  when predicting whether the observed temperature exceeds \code{th}.}
  \item{CondClim}{Probability forecast issued by a conditional climatological forecast distribution,
  which is different in each \code{state}, when predicting whether the observed temperature
  exceeds \code{th}.}
  \item{PP}{Probability forecast issued by a statistically post-processed forecast distribution
  when predicting whether the observed temperature exceeds \code{th}.}
  \item{CondPP}{Probability forecast issued by a conditional statistical post-processing method,
  which depends on \code{state}, when predicting whether the observed temperature
  exceeds \code{th}.}
  \item{state}{Integer from 1 to 5 corresponding to the state of the forecast system.}
  \item{th}{Temperature threshold from 5 to 30 degrees celcius.}
}
}
\usage{
data("noisy_data")
}
\description{
Simulated data designed to resemble probability forecasts issued by MeteoSwiss's
COSMO-E ensemble prediction system when predicting whether the temperature in
Switzerland will exceed certain thresholds.
}
\details{
Forecasts are available for 5 prediction systems: MeteoSwiss's COSMO-E ensemble prediction system,
a climatological forecaster, a state-dependent climatological forecaster, a statistically post-processed
forecaster, and a state-dependent statistically post-processed forecaster.

All forecasts are probabilities between 0 and 1. The forecasts have all been mapped to a finite number
of possible values 1/23, 2/23, ..., 22/23. The corresponding observations are either 0 or 1.

In each forecasting case, \code{state} represents the conditions surrounding this forecast.
These states could represent, for example, different seasons, weather regimes, or spatial regions.
In this case, there are five states, representing groupings of weather stations in Switzerland
depending on their altitude.

Each prediction system issues forecasts for 26 thresholds from 5 to 30 degrees celcius.
There are 38908 separate forecast cases, resulting in 38908*26=1011608 rows in the dataframe.
}
\keyword{datasets}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DiversityOccu.R
\name{model.diversity}
\alias{model.diversity}
\title{Find the best GLM model explaining the alpha divesity of the species}
\usage{
model.diversity(DivOcc)
}
\arguments{
\item{DivOcc}{is an object returned by the divesityoccu function of this package}
}
\value{
An object with the best fitted model, the coefficients of that model
and a table with the top 5 fitted models ranked by AICc
}
\description{
This function takes a diversityoccu object and heuristically searches for the
glm that best explains the alpha diversity of the modelled species.
}
\details{
This function fits every first order glm possible and ranks them by AICc
}
\examples{
data("BatOccu")
data("Dailycov")
data("sampling.cov")
x <-diversityoccu(pres = BatOccu, sitecov = sampling.cov, obscov = Dailycov, spp = 17, form = ~ Julian + Meanhum + Meantemp + sdhum + sdtemp ~  Burn.intensity.soil + I(Burn.intensity.soil^2) + Burn.intensity.Canopy + I(Burn.intensity.Canopy^2) + Burn.intensity.basal + I(Burn.intensity.basal^2))
model.diversity(x)
}
\author{
Derek Corcoran <derek.corcoran.barrios@gmail.com>
}
\seealso{
\code{\link[DiversityOccupancy]{diversityoccu}}
}

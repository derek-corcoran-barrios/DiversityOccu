#' summarizes the diversity and abundances predicted from a diversityoccupancy
#' object

#' @examples
#' data("BatOccu")
#' data("Dailycov")
#' data("sampling.cov")
#' x <-diversityoccu(pres = BatOccu, sitecov = sampling.cov, obscov = Dailycov,
#' spp = 17, form = ~ Julian + Meanhum + Meantemp + sdhum + sdtemp ~
#' Burn.intensity.soil + I(Burn.intensity.soil^2) + Burn.intensity.Canopy +
#' I(Burn.intensity.Canopy^2) + Burn.intensity.basal +
#' I(Burn.intensity.basal^2))
#'
#' summary(x)
#'
#' @export
#' @seealso \code{\link[DiversityOccupancy]{diversityoccu}}
#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>

summary.diversityoccupancy <- function (x) {
  summary(x$species)
}

#' summarizes the diversity and abundances predicted from a diversityoccupancy
#' object

#' @examples
#' data("BatOccu")
#' data("Dailycov")
#' data("sampling.cov")
#' x <-diversityoccu(pres = BatOccu, sitecov = sampling.cov, obscov = Dailycov,
#' spp = 17, form = ~ Julian + Meanhum + Meantemp + sdhum + sdtemp ~
#' Burn.intensity.soil + I(Burn.intensity.soil^2) + Burn.intensity.Canopy +
#' I(Burn.intensity.Canopy^2) + Burn.intensity.basal +
#' I(Burn.intensity.basal^2))
#' y <- model.diversity(x, method = "g", squared = TRUE)
#' summary(y)
#'
#' @export
#' @seealso \code{\link[DiversityOccupancy]{diversityoccu}}
#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>

summary.modeldiversity <- function (x) {
  x$Table
  }

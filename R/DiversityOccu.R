#' Calculates alpha diversity from multiple species occupancy data
#'
#' This function takes a data.frame with multiple presence absence-data from various species in different sites, covariates of each site to calculate occupancy, variables specific to sampling days to calculate probability of detection, and it calculates the alpha diversity for each site.
#'
#' @param pres a data.frame where rows are the sites and columns are a series of presence-absence observation from multiple species, every specie needs to have the same number of observations.
#' @param sitecov a data.frame where every row is a site, and every column is a measurement of that site, such as elevation or slope, this covariates are usually more constant.
#' @param obscov a list where every element is a data frame with the daily covariates for each site, that is a measurement for each day, such as average temperature of a day, this covariates are usually very .
#' @param spp the number of species in the pres data.frame
#' @param form a formula in the format ~ obscov ~ sitcov, the first arguments will be used to calculate probability of detection and the second part the occupancy.
#' @return A list with the fitted models for each specie and the calculated Alpha diversity for each site.
#' @details
#' This function fits the latent abundance mixture model described in Royle and Nichols (2003), to calculate the abundance of every species in each site, the using that abundance it calculates the alpha diversity index for each site based on that abundance.
#' @examples
#' data("BatOccu")
#' data("Dailycov")
#' data("sampling.cov")
#' diversityoccu(pres = BatOccu, sitecov = sampling.cov, obscov = Dailycov, spp = 17, form = ~ Julian + Meanhum + Meantemp + sdhum + sdtemp ~  Burn.intensity.soil + I(Burn.intensity.soil^2) + Burn.intensity.Canopy + I(Burn.intensity.Canopy^2) + Burn.intensity.basal + I(Burn.intensity.basal^2))
#' @seealso \code{\link[vegan]{diversity}}
#' @export
#' @importFrom vegan diversity
#' @importFrom unmarked occuRN
#' @importFrom unmarked unmarkedFrameOccu
#' @importMethodsFrom unmarked predict
#'

#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>

diversityoccu<- function(pres, sitecov, obscov, spp, form) {

  secuencia <- c(1:spp)*(ncol(pres)/spp)
  secuencia2<-secuencia-(secuencia[1]-1)

  models <- list()
  div <- list()

  for(i in 1:length(secuencia)) {
    models[[i]] <-c(secuencia2[i]:secuencia[i])
    models[[i]] <- pres[, models[[i]]]
    models[[i]] <- unmarkedFrameOccu(y = models[[i]], siteCovs = sitecov, obsCovs = obscov)
    models[[i]] <- occuRN(form, models[[i]])
    div[[i]] <- predict(models[[i]], type = "state", newdata = sitecov)$Predicted
    div<- as.data.frame(div)
    h<- diversity(div)
  }

  result <- list(Covs = sitecov, models = models, Diversity = h)
  return(result)
}

# You need the suggested package for this function
my_fun <- function(a, b) {
  if (!requireNamespace("pkg", quietly = TRUE)) {
    stop("Pkg needed for this function to work. Please install it.",
         call. = FALSE)
  }
}

#' Find the best GLM model explaining the alpha divesity of the species
#'
#' This function takes a diversityoccu object and heuristically searches for the
#' glm that best explains the alpha diversity of the modelled species.
#'
#' @param DivOcc is an object returned by the divesityoccu function of this package
#' @return An object with the best fitted model, the coefficients of that model
#' and a table with the top 5 fitted models ranked by AICc
#' @details
#' This function fits every first order glm possible and ranks them by AICc
#' @examples
#' data("BatOccu")
#' data("Dailycov")
#' data("sampling.cov")
#' x <-diversityoccu(pres = BatOccu, sitecov = sampling.cov, obscov = Dailycov, spp = 17, form = ~ Julian + Meanhum + Meantemp + sdhum + sdtemp ~  Burn.intensity.soil + I(Burn.intensity.soil^2) + Burn.intensity.Canopy + I(Burn.intensity.Canopy^2) + Burn.intensity.basal + I(Burn.intensity.basal^2))
#' model.diversity(x)
#' @seealso \code{\link[DiversityOccupancy]{diversityoccu}}
#' @export
#' @importFrom glmulti glmulti
#' @importFrom glmulti weightable
#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>


model.diversity <- function(DivOcc){
  A <- cbind(DivOcc$Diversity, DivOcc$Covs)
  colnames(A)[1]<-"Diversity"
  B <- paste(names(DivOcc$Covs), "+")
  B <- toString(B)
  B <- gsub(",", " ", B)
  B <- paste("Diversity ~", B, collapse = " ")
  B <- as.formula(substr(B, 1, nchar(B)-1))
  B <- glm(B, data = A)
  C <- glmulti(B, level = 1, crit = "aicc", confsetsize = 5, plotty = FALSE, method = "g")
  Best.model <- C@formulas[[1]]
  Table <- weightable(C)
  d<-summary(glm(Best.model, data = A))
  result <- list(Best_model = Best.model, Table = Table, coeff = d)
  return(result)
}

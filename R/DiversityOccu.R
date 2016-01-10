#' Calculates alpha diversity from multiple species occupancy data
#'
#' This function takes a data.frame with multiple presence absence-data from
#' various species in different sites, covariates of each site to calculate
#' occupancy, variables specific to sampling days to calculate probability of
#' detection, and it calculates the alpha diversity for each site.
#'
#' @param pres a data.frame where rows are the sites and columns are a series of
#' presence-absence observation from multiple species, every species needs to
#' have the same number of observations.
#' @param sitecov a data.frame where every row is a site, and every column is a
#' measurement of that site, such as elevation or slope, this covariates are
#' usually more constant.
#' @param obscov a list where every element is a data frame with the daily
#' covariates for each site, that is a measurement for each day, such as average
#' temperature of a day, this covariates are usually very .
#' @param spp the number of species in the pres data.frame
#' @param form a formula in the format ~ obscov ~ sitcov, the first arguments
#' will be used to calculate probability of detection and the second part the
#' occupancy.
#' @param index Diversity index, one of "shannon", "simpson" or "invsimpson".
#' @return A list with the fitted models for each species and the calculated
#' Alpha diversity for each site.
#' @details
#' This function fits the latent abundance mixture model described in Royle and
#' Nichols (2003), to calculate the abundance of every species in each site, the
#' using that abundance it calculates the alpha diversity index for each site
#' based on that abundance.
#' @examples
#' data("BatOccu")
#' data("Dailycov")
#' data("sampling.cov")
#' diversityoccu(pres = BatOccu, sitecov = sampling.cov, obscov = Dailycov,
#' spp = 17, form = ~ Julian + Meanhum + Meantemp + sdhum + sdtemp ~
#' Burn.intensity.soil + I(Burn.intensity.soil^2) + Burn.intensity.Canopy +
#' I(Burn.intensity.Canopy^2) + Burn.intensity.basal +
#' I(Burn.intensity.basal^2))
#' @seealso \code{\link[vegan]{diversity}}
#' @export
#' @importFrom vegan diversity
#' @importFrom unmarked occuRN
#' @importFrom unmarked unmarkedFrameOccu
#' @importMethodsFrom unmarked predict
#'

#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>

diversityoccu<- function(pres, sitecov, obscov, spp, form, index = "shannon") {

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
    h<- diversity(div, index)
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
#' @param DivOcc is an object returned by the divesityoccu function of this
#' package
#' @param method The method to be used to explore the candidate set of models.
#' If "h" an exhaustive screening is undertaken. If "g" the genetic algorithm is
#' employed (recommended for large candidate sets). If "l", a very fast
#' exhaustive branch-and-bound algorithm is used. Package leaps must then be
#' loaded, and this can only be applied to linear models with covariates and no
#' interactions.
#' @param confsetsize	The number of models to be looked for, i.e. the size of
#' the returned confidence set.
#' @return An object with the best fitted model, the coefficients of that model,
#' a table with the top 5 fitted models ranked by AICc and the data used for the
#' model
#' @details
#' This function fits every first order glm possible and ranks them by AICc
#' @examples
#' data("BatOccu")
#' data("Dailycov")
#' data("sampling.cov")
#' x <-diversityoccu(pres = BatOccu, sitecov = sampling.cov, obscov = Dailycov,
#' spp = 17, form = ~ Julian + Meanhum + Meantemp + sdhum + sdtemp ~
#' Burn.intensity.soil + I(Burn.intensity.soil^2) + Burn.intensity.Canopy +
#' I(Burn.intensity.Canopy^2) + Burn.intensity.basal +
#' I(Burn.intensity.basal^2))
#' y <- model.diversity(x, method = "g")
#' y$Table
#' y
#' @seealso \code{\link[DiversityOccupancy]{diversityoccu}}
#' @export
#' @importFrom glmulti glmulti
#' @importFrom glmulti weightable
#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>


model.diversity <- function(DivOcc, method = "h", confsetsize = 5){
  A <- cbind(DivOcc$Diversity, DivOcc$Covs)
  colnames(A)[1]<-"Diversity"
  B <- paste(names(DivOcc$Covs), "+")
  B <- toString(B)
  B <- gsub(",", " ", B)
  B <- paste("Diversity ~", B, collapse = " ")
  B <- as.formula(substr(B, 1, nchar(B)-1))
  B <- glm(B, data = A)
  C <- glmulti(B, level = 1, crit = "aicc", confsetsize = confsetsize, plotty = FALSE, method = method)
  Best.model <- C@formulas[[1]]
  Table <- weightable(C)
  Table$Delta.AICc <- Table[,2]-Table[1,2]
  d<-summary(glm(Best.model, data = A))
  result <- list(Best_model = Best.model, Table = Table, coeff = d, dataset= A)
  return(result)
}

# You need the suggested package for this function
my_fun <- function(a, b) {
  if (!requireNamespace("pkg", quietly = TRUE)) {
    stop("Pkg needed for this function to work. Please install it.",
         call. = FALSE)
  }
}

#' plot the response of the calculated alpha diversity to the change of a
#' particular variable
#'
#' This function takes a model.diversity object and one of the variables used to
#' predict the
#' @examples
#' data("BatOccu")
#' data("Dailycov")
#' data("sampling.cov")
#' x <-diversityoccu(pres = BatOccu, sitecov = sampling.cov, obscov = Dailycov,
#' spp = 17, form = ~ Julian + Meanhum + Meantemp + sdhum + sdtemp ~
#' Burn.intensity.soil + I(Burn.intensity.soil^2) + Burn.intensity.Canopy +
#' I(Burn.intensity.Canopy^2) + Burn.intensity.basal +
#' I(Burn.intensity.basal^2))
#' y <- model.diversity(x, method = "g")
#' response.plot(y, Burn.intensity.soil)
#' response.plot(y, Existing.vegetation)
#' @export
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_ribbon
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_line
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 labs
#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>

response.plot<- function(model, variable){
  a<-data.frame(matrix(rep(colMeans(model$dataset), each=length(model$dataset[,1])), nrow = length(model$dataset[,1]), ncol = ncol(model$dataset)))
  colnames(a)<-colnames(model$dataset)
  maxval<-apply(model$dataset,2,max)
  minval<-apply(model$dataset,2,min)
  newdata<- seq(from = minval[colnames(a)== as.character(substitute(variable))], to = maxval[colnames(a)== as.character(substitute(variable))], along.with = model$dataset[,1])
  a[colnames(a)== as.character(substitute(variable))] <- newdata
  b<-predict(glm(model$Best_model, data= model$dataset), newdata = a, se.fit = TRUE)
  c<- data.frame(preditction = b$fit, upper = (b$fit + b$se), lower = (b$fit - b$se), dependent = a[colnames(a)== as.character(substitute(variable))])
  result <- ggplot(c, aes(x= c[,4], y = c[,1])) + geom_ribbon(aes(ymax= c[,2], ymin = c[,3]), fill = "grey") + geom_line() + theme_bw() + theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank()) + labs(x = as.character(substitute(variable)), y = "Diversity")
  return(result)
}

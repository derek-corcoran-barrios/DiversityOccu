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
#' @param dredge default = FALSE, if TRUE, for each species, the best occupancy
#' model will be determined by fitting all possible models and rankin by AICc.
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
#' @seealso \code{\link[DiversityOccupancy]{model.diversity}}
#' @seealso \code{\link[DiversityOccupancy]{response.plot}}
#' @export
#' @importFrom vegan diversity
#' @importFrom unmarked occuRN
#' @importFrom unmarked unmarkedFrameOccu
#' @importMethodsFrom unmarked predict
#' @importFrom MuMIn dredge
#' @importFrom MuMIn get.models
#' @importFrom MuMIn AICc

#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>
#' @author Nicole L. Michel

diversityoccu<- function(pres, sitecov, obscov, spp, form, index = "shannon", dredge = FALSE) {

  secuencia <- c(1:spp)*(ncol(pres)/spp)
  secuencia2<-secuencia-(secuencia[1]-1)

  models <- list()
  data<- list()
  div <- list()
  dredged <- list ()
  if (dredge == FALSE){

    for(i in 1:length(secuencia)) {
      data[[i]] <-c(secuencia2[i]:secuencia[i])
      data[[i]] <- pres[, data[[i]]]
      models[[i]] <- unmarkedFrameOccu(y = data[[i]], siteCovs = sitecov, obsCovs = obscov)
      models[[i]] <- occuRN(form, models[[i]])
      div[[i]] <- predict(models[[i]], type = "state", newdata = sitecov)$Predicted
      div<- as.data.frame(div)
      h<- diversity(div, index)
    }
  }

  else if (dredge==TRUE) {
    for(i in 1:length(secuencia)) {
      require(MuMIn)
      require(unmarked)
      data[[i]] <-c(secuencia2[i]:secuencia[i])
      data[[i]] <- pres[, data[[i]]]
      #data is a list of class unmarkedFrames from package unmarked.
      # NM: write to the global environment so he data won't be "lost"
      data2 <<- unmarkedFrameOccu(y = data[[i]], siteCovs = sitecov, obsCovs = obscov)
      #uses the data list above to fit one model per specie
      models[[i]] <- occuRN(form, data2)
      #selects models
      # NM: saved this to dredged object rather than overwriting models object
      dredged[[i]] <- dredge(models[[i]], data2)
      #select the first model
      models[[i]] <- get.models(dredged[[i]], 1)[[1]]
      #predictions for the best model
      div[[i]] <- predict(models[[i]], type = "state", newdata = sitecov)$Predicted
      div<- as.data.frame(div)
      colnames(div) = paste("species",c(1:ncol(div)), sep =".")
      data[[i]] <- data2
      h<- diversity(div, index)
    }
    # remove temporary data file from the global environment
    rm(data2, pos=".GlobalEnv")
  }

  result <- list(Covs = sitecov, models = models, Diversity = h)
  return(result)
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
#' @param delta	The number of models that will be returned will be the ones that
#' have a maximum AICc difference with the top model equal to delta.
#' @param squared, if FALSE (Default), only GLMs with linear components will be
#' evaluated; If TRUE, GLMs with both linear and quadratic components will be evaluated.
#' WARNING if squared is TRUE, the number of parameters duplicates and the models
#' grow exponentially, this may result in to many variables for a CPU to compute.
#' @return An object with the best fitted model, the coefficients of that model,
#' a table with the top 5 fitted models ranked by AICc and the data used for the
#' model
#' @details
#' This function fits every first order glm possible and ranks them by AICc.
#' @examples
#' #To fit and explore the only the linear components of the model
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
#'
#' #To add the quadratic components of models
#'
#' data("BatOccu")
#' data("Dailycov")
#' data("sampling.cov")
#' x <-diversityoccu(pres = BatOccu, sitecov = sampling.cov[,1:8], obscov = Dailycov,
#' spp = 17, form = ~ Julian + Meanhum + Meantemp + sdhum + sdtemp ~
#' Burn.intensity.soil + I(Burn.intensity.soil^2) + Burn.intensity.Canopy +
#' I(Burn.intensity.Canopy^2) + Burn.intensity.basal +
#' I(Burn.intensity.basal^2))
#' y <- model.diversity(x, method = "g", squared = TRUE)
#' y$Table
#' y
#'
#' @seealso \code{\link[DiversityOccupancy]{diversityoccu}}
#' @seealso \code{\link[DiversityOccupancy]{response.plot}}
#' @export
#' @importFrom glmulti glmulti
#' @importFrom glmulti weightable
#' @importFrom dplyr filter
#' @importFrom qpcR akaike.weights
#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>


model.diversity <- function(DivOcc, method = "h", delta = 2, squared = FALSE){
  A <- cbind(DivOcc$Diversity, DivOcc$Covs)
  colnames(A)[1]<-"Diversity"
  B <- paste(names(DivOcc$Covs), "+")
  B <- toString(B)
  B <- gsub(",", " ", B)
  if (squared == TRUE) {
    C <- paste("I(", names(DivOcc$Covs), "^2) +")
    C <- toString(C)
    C <- gsub (",", " ", C)
    B <- paste("Diversity ~", B, C, collapse = " ")
  }
  else if (squared == FALSE) {
    B <- paste("Diversity ~", B, collapse = " ")
  }

  B <- as.formula(substr(B, 1, nchar(B)-1))
  B <- glm(B, data = A)
  D <- glmulti(B, level = 1, crit = "aicc", plotty = FALSE, method = method)
  Best.model <- D@formulas[[1]]
  Table <- weightable(D)
  Table$Delta.AICc <- Table[,2]-Table[1,2]
  Table <- filter(Table, Delta.AICc < delta)
  Table$weights <- akaike.weights(Table$aicc)$weights
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
#' predict the alpha diversity index, and makes a plot showing the response of
#' the diversity index against the selected variable. This function automatically
#' limits the values of that variable to the maximum and minimum values of the
#' dataset.
#' @param model A result from the model.diversity function.
#' @param variable The variable of which the response is to be ploted.
#' @return a ggplot object plotting the alpha diversity response to the selected
#' variable.
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
#' @seealso \code{\link[DiversityOccupancy]{diversityoccu}}
#' @seealso \code{\link[DiversityOccupancy]{model.diversity}}
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
  A<-data.frame(matrix(rep(colMeans(model$dataset), each=length(model$dataset[,1])), nrow = length(model$dataset[,1]), ncol = ncol(model$dataset)))
  colnames(A)<-colnames(model$dataset)
  maxval<-apply(model$dataset,2,max)
  minval<-apply(model$dataset,2,min)
  newdata<- seq(from = minval[colnames(A)== as.character(substitute(variable))], to = maxval[colnames(A)== as.character(substitute(variable))], along.with = model$dataset[,1])
  A[colnames(A)== as.character(substitute(variable))] <- newdata
  B<-predict(glm(model$Best_model, data= model$dataset), newdata = A, se.fit = TRUE)
  C<- data.frame(preditction = B$fit, upper = (B$fit + B$se), lower = (B$fit - B$se), dependent = A[colnames(A)== as.character(substitute(variable))])
  result <- ggplot(C, aes(x= C[,4], y = C[,1])) + geom_ribbon(aes(ymax= C[,2], ymin = C[,3]), fill = "grey") + geom_line() + theme_bw() + theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank()) + labs(x = as.character(substitute(variable)), y = "Diversity")
  return(result)
}

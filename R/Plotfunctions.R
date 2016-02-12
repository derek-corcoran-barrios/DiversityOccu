
#' plot the response of the calculated alpha diversity to the change of a
#' particular variable
#'
#' This function takes a modeldiversity object and one of the variables used to
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
#' BatDiversity <-diversityoccu(pres = BatOccu, sitecov = sampling.cov[,1:8],
#' obscov = Dailycov,spp = 17, form = ~ Julian + Meanhum + Meantemp + sdhum +
#' sdtemp ~ Burn.intensity.soil + I(Burn.intensity.soil^2) +
#' Burn.intensity.Canopy + I(Burn.intensity.Canopy^2) + Burn.intensity.basal +
#' I(Burn.intensity.basal^2))
#'
#' #Select the best model that explains diversity using genetic algorithms
#' set.seed(123)
#' glm.Batdiversity <- model.diversity(x, method = "g")
#'
#' #see the best models
#'
#' summary(glm.Batdiversity)
#'
#' #plot the response of diversity to individual variables
#'
#' plot(glm.Batdiversity, Burn.intensity.soil)
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

plot.modeldiversity<- function(model, variable){
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

#' plots a histogram of diversity of a predicted diversityoccupancy object

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
#' plot (x)
#'
#' @export
#' @seealso \code{\link[DiversityOccupancy]{diversityoccu}}
#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>

plot.diversityoccupancy <- function (x) {
  hist(x$Diversity, xlab = "Alpha diversity", ylab = "Number of sites", main = NULL)
}

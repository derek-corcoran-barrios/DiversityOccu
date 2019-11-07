#' plot the response of an occupancy model to the change of aparticular variable
#'
#' This function takes a batchoccupancy object and one of the variables used to
#' predict occupacy, and makes a plot showing the response of occupancyt against
#' the selected variable. This function automatically limits the values of that
#' variable to the maximum and minimum values of the dataset.
#' @param batch A result from the batchoccu function.
#' @param spp The species number of which response is going to be ploted.
#' @param variable The variable of which the response is to be ploted.
#' @param N The number of points used to make the plot (default is 20) the greater
#' the number, the clearer will be the response, but it will make the plot slower
#' @return a ggplot object plotting the alpha diversity response to the selected
#' variable.
#' @examples
#' \dontrun{
#' data("IslandBirds")
#' data("Daily_Cov")
#' data("siteCov")
#' BirdOccupancy <-batchoccu(pres = IslandBirds, sitecov = siteCov, obscov =
#' Daily_Cov, spp = 5, form = ~ Day + Wind + Noise + Clouds ~
#' Elev + AgroFo + SecVec + Wetland)
#' #plot the response of occupancy to individual variables for species 4 and 5
#'
#' responseplot.occu(batch = BirdOccupancy, spp = 4, variable = "Elev")
#'
#'
#' responseplot.occu(batch = BirdOccupancy, spp = 5, variable = "Elev")
#' }
#' @seealso \code{\link[DiversityOccupancy]{batchoccu}}
#' @importFrom dplyr "%>%"
#' @importFrom dplyr bind_cols
#' @importFrom dplyr one_of
#' @importFrom dplyr pull
#' @importFrom dplyr quo
#' @importFrom dplyr select
#' @importFrom dplyr select_if
#' @importFrom dplyr summarize_all
#' @importFrom dplyr syms
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_errorbar
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_ribbon
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 sym
#' @importFrom ggplot2 theme_classic
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 ylim
#' @importFrom purrr map
#' @importFrom stringr str_count
#' @importFrom stringr str_remove_all
#' @importFrom stringr str_split
#' @importFrom stringr str_trim
#' @importFrom rlang ":="
#' @importFrom rlang .data
#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>
#' @export

responseplot.occu <- function (batch, spp, variable, N = 50){
  upper <- lower <- Predicted <- NULL

  Vars <- batch$models[[spp]]@formula %>% str_split(pattern = "~", simplify = T)
  Vars <- Vars[3] %>% str_split(pattern = "\\+", simplify = T) %>% as.character() %>% str_trim() %>% str_remove_all("1")
  Vars <- Vars[str_count(Vars) != 0]

  model_data <- batch$Covs %>% select(Vars)

  stopifnot(variable %in% Vars)

  all_vars <- model_data %>% dplyr::select(-one_of(variable))
  num_vars <- all_vars %>% select_if(is.numeric) %>% summarize_all(mean)
  cat_vars <- all_vars %>% select_if(Negate(is.numeric)) %>% purrr::map(unique)

  resp_var <- model_data %>% pull(variable)
  if(is.numeric(resp_var)) {
    resp_vals <- seq(min(resp_var), max(resp_var), length.out=N)
  } else {
    resp_vals <- unique(resp_var)
  }

  new_data <- tidyr::crossing(num_vars, !!!cat_vars, !!variable:=resp_vals)
  if(nrow(new_data) < 1){
    new_data <- data.frame(resp_vals)
    colnames(new_data) <- variable
  }

  pred <- predict(batch$models[[spp]], newdata = as.data.frame(new_data), se_fit=TRUE, type= "state") %>% bind_cols(new_data)

  ## Plot the response
  my_aes <- aes(x= !!sym(variable), y = Predicted)
  line_aes <- aes(x= !!sym(variable), y = Predicted)
  if (length(cat_vars)==1) {
    my_aes[["fill"]] <- sym(names(cat_vars))
    line_aes[["colour"]] <- sym(names(cat_vars))
  } else if (length(cat_vars)>1) {
    my_aes[["fill"]] <- quo(interaction(!!!syms(names(cat_vars))))
    line_aes[["colour"]] <- quo(interaction(!!!syms(names(cat_vars))))
  }
  range_aes <- aes(ymax= upper, ymin = lower)
  result <- ggplot(pred, my_aes) + theme_classic() + ylab("Occupancy") + ylim(c(0,1))
  if(is.numeric(resp_var)) {
    (if (length(cat_vars)>0) {
      result <- result + geom_ribbon(range_aes, alpha = 0.5) + geom_line(line_aes)
    } else {
      result <- result + geom_ribbon(range_aes, fill="grey", alpha = 0.5) + geom_line()
    })
  } else {
    result +
      geom_errorbar(range_aes) +
      geom_point()
  }

}

#' plot the response of an abundance model to the change of aparticular variable
#'
#' This function takes a diversityoccupancy object and one of the variables used
#' to predict abundance, and makes a plot showing the response of occupancyt
#' against the selected variable. This function automatically limits the values
#' of that variable to the maximum and minimum values of the dataset.
#' @param batch A result from the diversityoccu function.
#' @param spp The species number of which response is going to be ploted.
#' @param variable The variable of which the response is to be ploted.
#' @return a ggplot object plotting the alpha diversity response to the selected
#' variable.
#' @examples
#' \dontrun{
#' data("IslandBirds")
#' data("Daily_Cov")
#' data("siteCov")
#'
#' #Model the abundance for 5 bird species and calculate alpha diversity from that
#'
#' BirdDiversity <-diversityoccu(pres = IslandBirds, sitecov = siteCov,
#' obscov = Daily_Cov,spp = 5, form = ~ Day + Wind + Time + Rain +
#' Noise ~ Elev + AgroFo + SecVec + Wetland + Upland)
#'
#' #plot the response of abundance to individual variables for species 4, 11
#'
#' responseplot.abund(batch = BirdDiversity, spp = 4, variable = Elev)
#'
#' responseplot.abund(batch = BirdDiversity, spp = 11, variable = Elev)
#' }
#' @export
#' @seealso \code{\link[DiversityOccupancy]{batchoccu}}
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_ribbon
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_line
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 ylim
#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>

responseplot.abund <- function(batch, spp, variable){
  upper <- lower <- NULL # Setting the variables to NULL first
  A<-data.frame(matrix(rep(colMeans(batch$Covs), each=length(batch$Covs[,1])), nrow = length(batch$Covs[,1]), ncol = ncol(batch$Covs)))
  colnames(A)<-colnames(batch$Covs)
  maxval<-apply(batch$Covs,2,max)
  minval<-apply(batch$Covs,2,min)
  newdata<- seq(from = minval[colnames(A)== as.character(substitute(variable))], to = maxval[colnames(A)== as.character(substitute(variable))], along.with = batch$Covs[,1])
  A[colnames(A)== as.character(substitute(variable))] <- newdata
  B<-predict(batch$models[[spp]], type = "state", newdata = A)
  DF<- data.frame(preditction = B$Predicted, upper = (B$Predicted + B$SE), lower = (B$Predicted - B$SE), dependent = A[colnames(A)== as.character(substitute(variable))])
  result <- ggplot(DF, aes(x= DF[,4], y = DF[,1])) + geom_ribbon(aes(ymax= upper, ymin = lower), fill = "grey") + geom_line() + theme_bw() + theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank()) + labs(x = as.character(substitute(variable)), y = "Abundance") + ylim(c(min(c(DF$lower, 0)),max(DF$upper)))
  return(result)
}

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
#' \dontrun{
#' #Load the data
#' data("IslandBirds")
#' data("Daily_Cov")
#' data("siteCov")
#'
#' #Model the abundance for 5 bird species and calculate alpha diversity from that
#'
#' BirdDiversity <-diversityoccu(pres = IslandBirds, sitecov = siteCov,
#' obscov = Daily_Cov,spp = 5, form = ~ Day + Wind + Time + Rain +
#' Noise ~ Elev + AgroFo + SecVec + Wetland + Upland)
#'
#' #Select the best model that explains diversity using genetic algorithms
#' set.seed(123)
#' glm.Birdiversity <- model.diversity(BirdDiversity, method = "g")
#'
#' #see the best models
#'
#' glm.Birdiversity$Best.model
#'
#' #plot the response of diversity to individual variables
#'
#' plot(glm.Birdiversity, elev)
#' }
#' @export
#' @seealso \code{\link[DiversityOccupancy]{diversityoccu}}
#' @seealso \code{\link[DiversityOccupancy]{model.diversity}}
#' @importFrom stats glm
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

responseplot.diver<- function(model, variable){
  A<-data.frame(matrix(rep(colMeans(model$dataset), each=length(model$dataset[,1])), nrow = length(model$dataset[,1]), ncol = ncol(model$dataset)))
  colnames(A)<-colnames(model$dataset)
  maxval<-apply(model$dataset,2,max)
  minval<-apply(model$dataset,2,min)
  newdata<- seq(from = minval[colnames(A)== as.character(substitute(variable))], to = maxval[colnames(A)== as.character(substitute(variable))], along.with = model$dataset[,1])
  A[colnames(A)== as.character(substitute(variable))] <- newdata
  B<-predict(glm(model$Best_model, data= model$dataset), newdata = A, se.fit = TRUE)
  DF<- data.frame(preditction = B$fit, upper = (B$fit + B$se), lower = (B$fit - B$se), dependent = A[colnames(A)== as.character(substitute(variable))])
  result <- ggplot(DF, aes(x= DF[,4], y = DF[,1])) + geom_ribbon(aes(ymax= DF[,2], ymin = DF[,3]), fill = "grey") + geom_line() + theme_bw() + theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank()) + labs(x = as.character(substitute(variable)), y = "Diversity")
  return(result)
}

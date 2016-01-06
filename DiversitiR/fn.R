ncol(BatOccu)
spp<-6
form <- ~ Julian + Meanhum + Meantemp + sdhum + sdtemp ~  Burn.intensity.soil + I(Burn.intensity.soil^2) + Burn.intensity.Canopy + I(Burn.intensity.Canopy^2) + Burn.intensity.basal + I(Burn.intensity.basal^2)
ncol(BatOccu)/spp

secuencia <- c(1:spp)*(ncol(BatOccu)/spp)
secuencia2<-secuencia-(secuencia[1]-1)

models <- list()
div <- list()

for(i in 1:length(secuencia))
{
  models[[i]] <-c(secuencia2[i]:secuencia[i])
  models[[i]] <- BatOccu[,models[[i]]]
  models[[i]] <- unmarkedFrameOccu(y = models[[i]], siteCovs =sampling.cov2, obsCovs=Dailycov3)
  models[[i]] <- occuRN(form, models[[i]])
  div[[i]] <- predict(models[[i]], type = "state", newdata=sampling.cov2)$Predicted
  div<- as.data.frame(div)
  h<- diversity(div)
}

result <- list(models=models,divers=h)


diversityoccu<- function(pres, sitecov, obscov, spp, form) {
  
  secuencia <- c(1:spp)*(ncol(pres)/spp)
  secuencia2<-secuencia-(secuencia[1]-1)
  
  models <- list()
  div <- list()
  
  for(i in 1:length(secuencia))
  {
    models[[i]] <-c(secuencia2[i]:secuencia[i])
    models[[i]] <- pres[,models[[i]]]
    models[[i]] <- unmarkedFrameOccu(y = models[[i]], siteCovs =sitecov, obsCovs=obscov)
    models[[i]] <- occuRN(form, models[[i]])
    div[[i]] <- predict(models[[i]], type = "state", newdata=sitecov)$Predicted
    div<- as.data.frame(div)
    h<- diversity(div)
  }
  
  result <- list(models=models,Diversity=h)
  return(result)
}

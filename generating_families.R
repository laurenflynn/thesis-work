set.seed(777)

# Load Libraries ----

setwd("../PedUtils/R")
files.sources = list.files()
sapply(files.sources, source) #loads in all the PedUtils functions
setwd("../../thesis-work")
source("family_filters.R")
library(plyr) #need to load plyr before dplyr
library(dplyr)
library(truncnorm)
library(PanelPRO)
library(tidyverse)
library(MASS)
library(parallel)
library(foreach)
library(doParallel)

generateFamilies <- function(mlh1Freq, numberFamilies,cores){
  # Set MLH1 Frequency to user provided val ----
  load("~/PanelPRO/data/PanelPRODatabase.rda")
  PanelPRODatabase$AlleleFrequency[["MLH1_anyPV",1]] <- mlh1Freq
  PanelPRODatabase$AlleleFrequency[["MLH1_anyPV",2]] <- mlh1Freq
  PanelPRODatabase$AlleleFrequency[["MLH1_anyPV",3]] <- mlh1Freq
  save(PanelPRODatabase, file = "~/PanelPRO/data/PanelPRODatabase.rda")
  
  
  # Generate Families ----
  families = list()
  probandIDS = c()
  probandMLH1Status = c()

  cl <- makeCluster(cores)
  registerDoParallel(cl)
  families <- foreach(i=1:numberFamilies, .packages =c("truncnorm","tidyverse")) %dopar%{
    setwd("../PedUtils/R")
    files.sources = list.files()
    sapply(files.sources, source) #loads in all the PedUtils functions
    setwd("../../thesis-work")
    # Cancers
    cancers = "Colorectal"
    # Genes
    genes = "MLH1"
    #family members
    # Paternal aunts, paternal uncles
    nSibsPatern =floor(rtruncnorm(n=2, mean=3, 3))
    # Maternal aunts, maternal uncles
    nSibsMatern = floor(rtruncnorm(n=2, mean=3, 3))
    # Sisters and brothers
    nSibs = floor(rtruncnorm(n=2, mean=3, 3))
    # We make the assumption that the number of sons and daughters for the
    # proband and all siblings, is the same. Nieces and nephews of the proband
    # are not sampled separately.
    nGrandchild = floor(rtruncnorm(n=2, mean=6, 2))
    nChild = floor(rtruncnorm(n=2, mean=3, 2))
    
    # Simulate family using `PedUtils` code
    fam = sim.runSimFam(nSibsPatern, nSibsMatern, nSibs, nChild,
                        PanelPRODatabase, genes, cancers,
                        includeGeno = TRUE, includeBiomarkers = TRUE)
    famDF = as.data.frame(fam)
    for(j in 1:nrow(famDF)){
      if(famDF[j,]$isAffAny == 0){
        famDF[j,]$AgeAny = NA    
      }
      if(famDF[j,]$isAffCOL == 0){
        famDF[j,]$AgeCOL = NA    
      }
    }
    proband = famDF %>% filter(isProband==1)
    probandIDS = c(probandIDS, proband$ID)
    probandMLH1Status = c(probandMLH1Status, proband$MLH1)
    # PanelPRO can be run on the simulated family
    #out = PanelPRO:::PanelPRO11(fam)
    families[[i]] = famDF
    #outputs = c(outputs, out)
  }
  
  probandIDS <- foreach(i=1:numberFamilies, .packages="tidyverse") %dopar%{
    proband <- families[[i]] %>% filter(isProband==1)
    probandID <- proband$ID
  }
  
  probandMLH1Status <- foreach(i=1:numberFamilies, .packages = "tidyverse") %dopar%{
    proband <- families[[i]] %>% filter(isProband==1)
    probandMLH1Status <- proband$MLH1
  }
  
  mlh1StatusFamilies <- families
  families <- foreach(i=1:numberFamilies, .export='removeProbandStatus') %dopar%{
    families[[i]]=removeProbandStatus(families[[i]])
  }
  
  firstDegree <- foreach(i=1:numberFamilies, .export = 'firstDegreeFamilyMembers') %dopar%{
    firstDegreeFamilyMembers(families[[i]])
  }
  
  maskedFamilies <- foreach(i=1:numberFamilies, .export = 'maskedInfoAmbry') %dopar%{
    maskedInfoAmbry(families[[i]])
  }
  # clusterExport(cl, c("families")) #not sure if we need this if we generate families in the cluster
  # families <- mclapply(families, removeProbandStatus, mc.cores = cores, mc.preschedule=FALSE) #removes mlh1 status from probands
  # firstDegree <- mclapply(families, firstDegreeFamilyMembers, mc.cores = cores, mc.preschedule=FALSE) #filters to first degree
  # maskedFamilies <- mclapply(families, maskedInfoAmbry, mc.cores = cores, mc.preschedule = FALSE)
  # stopCluster(cl)
  
  print("Full Families from PedUtils with MLH1 Status for Probands")
  describeFamilies(mlh1StatusFamilies)
  print("*******")
  print("Full Families from PedUtils with Unaffected Information Masked")
  describeFamilies(maskedFamilies)
  print("*******")
  print("First Degree Families")
  describeFamilies(firstDegree)
  
  # # cl <- makeCluster(cores)
  # # registerDoParallel(cl)
  # 
  # 
  outputsFull <- list ()
  fullCarrierRisk <- list()
  for(i in 1:length(families)){
    print(i)
    out = PanelPRO::PanelPRO(families[[i]], genes="MLH1", cancers="Colorectal")
    id = as.character(probandIDS[i])
    fullCarrierRisk[[i]] <- out$posterior.prob[[id]]$estimate[2]
    outputsFull[[i]] = out#DF
  }
  print("Ran PanelPRO on Full Families")
  # 
  # maskedFamiliesCarrierRisk <- foreach(i=1:numberFamilies, .packages="PanelPRO") %dopar% {
  #   out = PanelPRO::PanelPRO(maskedFamilies[[i]], genes="MLH1", cancers="Colorectal")
  #   id = as.character(probandIDS[i])
  #   out$posterior.prob[[id]]$estimate[2]
  # }
  outputsMasked <- list ()
  maskedCarrierRisk <- list()
  for(i in 1:length(families)){
    print(i)
    out = PanelPRO::PanelPRO(maskedFamilies[[i]], genes="MLH1", cancers="Colorectal")
    id = as.character(probandIDS[i])
    maskedCarrierRisk[[i]] <- out$posterior.prob[[id]]$estimate[2]
    outputsMasked[[i]] = out#DF
  }
  print("Ran PanelPRO with unaffected family members info masked")

  outputsFirstDegree <- list ()
  firstDegreeCarrierRisk <- list()
  for(i in 1:length(families)){
    print(i)
    out = PanelPRO::PanelPRO(firstDegree[[i]], genes="MLH1", cancers="Colorectal")
    id = as.character(probandIDS[i])
    firstDegreeCarrierRisk[[i]] <- out$posterior.prob[[id]]$estimate[2]
    outputsFirstDegree[[i]] = out#DF
  }
  print("Ran PanelPRO on first degree families")

  mlh1ProbandsTotal <- c()
  probandMLH1 <- list()
  for(i in 1:length(families)){
    if(nrow(mlh1StatusFamilies[[i]] %>% filter(MLH1 ==1) %>% filter(isProband==1) > 0)){
      mlh1ProbandsTotal <- c(mlh1ProbandsTotal, i)
      probandMLH1[[i]] = 1
    } 
    else{
      probandMLH1[[i]] = 0
    }}
  stopCluster(cl)


  summaryTableMLH1 <- cbind(format(unlist(fullCarrierRisk), scientific=FALSE), unlist(maskedFamiliesCarrierRisk), unlist(firstDegreeCarrierRisk), unlist(probandMLH1))
  summaryTableMLH1 <- as.data.frame(summaryTableMLH1)
  summaryTable <- cbind(1:numberFamilies, summaryTableMLH1)
  print(head(summaryTable))
  names(summaryTable) <- c("famID","fullCarrierRisk", "carrierRiskUnaffectedInfoMasked", "firstDegreeCarrierRisk", "probandMLH1Status")
  summaryTable$fullCarrierRisk <- as.numeric(summaryTable$fullCarrierRisk)
  summaryTable$carrierRiskUnaffectedInfoMasked <- as.numeric(summaryTable$carrierRiskUnaffectedInfoMasked)
  summaryTable$firstDegreeCarrierRisk <- as.numeric(summaryTable$firstDegreeCarrierRisk)

  visual = ggplot(data= summaryTable) +  geom_point(aes(x=famID, y=fullCarrierRisk, color="Full Carrier Risk")) +
    geom_point(aes(x=famNumber, y= carrierRiskUnaffectedInfoMasked, color = "Carrier Risk with Unaffected Info Masked")) +
    geom_point(aes(x=famNumber, y= firstDegreeCarrierRisk, color = "First Degree Carrier Risk")) +
    scale_y_log10() +
    labs(title= "MLH1 Carrier Risk", x= "Family ID", y="Carrier Risk")

  print(visual)

  fitMedianPolished <- rlm(fullCarrierRisk ~ carrierRiskUnaffectedInfoMasked, data = summaryTable)
  print(summary(fitMedianPolished))

  # finally, plot the data and the fitted curve
  medPolVis = ggplot(data = summaryTable, aes(x= carrierRiskUnaffectedInfoMasked, y= fullCarrierRisk)) +
    geom_point() +
    geom_line(aes(y = predict(fitMedianPolished)))+
    labs(x = "Carrier Risk Unaffected Info Masked", y= "Carrier Risk for Full Family Information", title = "Comparing Carrier Risk for Masked Info and Full Family Information with Loess Regression")

  print(medPolVis)

  fitLoess <- loess(fullCarrierRisk ~ carrierRiskUnaffectedInfoMasked, data = summaryTable, span = 0.2)
  summary(fitLoess) #RSE is 0.00397 (seems good but I'm not sure)

  # finally, plot the data and the fitted curve
  loessVis = ggplot(data = summaryTable, aes(x= carrierRiskUnaffectedInfoMasked, y= fullCarrierRisk)) +
    geom_point() +
    geom_line(aes(y = predict(fitLoess)))+
    labs(x = "Carrier Risk Unaffected Info Masked", y= "Carrier Risk for Full Family Information", title = "Comparing Carrier Risk for Masked Info and Full Family Information with Loess Regression")

  print(loessVis)
}

generateFamilies(0.2, 5, 5)

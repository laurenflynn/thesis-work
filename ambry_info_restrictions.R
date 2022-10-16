set.seed(777)
setwd("/Users/flynnlauren7/Desktop/Internships Summer 2022/Bayes Mendel/PedUtils/R")
files.sources = list.files()
sapply(files.sources, source)
library(dplyr)


families = list()
#outputs = c()
for (i in 1:10) {
  # Cancers
  cancers = "Colorectal"
  # Genes

  genes = "MLH1"

  # Paternal aunts, paternal uncles
  nSibsPatern =floor(rtruncnorm(n=2, mean=3, 2))
  # Maternal aunts, maternal uncles
  nSibsMatern = floor(rtruncnorm(n=2, mean=3, 2))
  # Sisters and brothers
  nSibs = floor(rtruncnorm(n=2, mean=3, 2))
  # We make the assumption that the number of sons and daughters for the
  # proband and all siblings, is the same. Nieces and nephews of the proband
  # are not sampled separately.
  nGrandchild = floor(rtruncnorm(n=2, mean=6, 2))
  nChild = floor(rtruncnorm(n=2, mean=3, 2))

  # Simulate family using `PedUtils` code
  fam = sim.runSimFam(nSibsPatern, nSibsMatern, nSibs, nChild,
                      PanelPRODatabase, genes, cancers,
                      includeGeno = FALSE, includeBiomarkers = TRUE)
  famDF = as.data.frame(fam)
  # PanelPRO can be run on the simulated family
  #out = PanelPRO:::PanelPRO11(fam)
  families[[i]] = famDF
  #outputs = c(outputs, out)

}



## Function to filter down to Ambry info ##
#### Proband info: sex, age, isDead==no, NA for riskMod and for interAge ###
ambryProbandInfo <- function(ped){
  proband = ped %>% filter(isProband==1)
  proband$riskmod <- NA
  proband$interAge <- NA #should these be character(0)? what does that mean?
  proband$isDead <- 0
  return(proband)
}

### relative info ###
#for simplicity, we will just say that sex can be inferred
#restricting to first degree family
#here genes is just MLH1
ambryRelativeInfo <- function(ped){
  #firstDegreeFam = firstDegreeFamilyMembers(ped)
  relatives = ped %>% filter(isProband==0)
  #print(relatives)
  relatives$curAge <- NULL
  relatives$isDead <- NULL
  relatives$riskMod <- NULL
  relatives$interAge <- NULL
  return(relatives)
}


firstDegreeFamilyMembers <- function(ped){
  proband <- ped %>% filter(isProband==1)
  firstDegree = proband
  if(!is.null(proband)){
    for(i in  1:nrow(ped)){
      relative=ped[i,]
      #siblings
      if(!is.null(relative$MotherID) && !is.null(proband$MotherID) && !is.null(relative$FatherID) && !is.null(proband$FatherID) && (relative$MotherID==proband$MotherID) && (relative$FatherID==proband$FatherID)){
        firstDegree <- rbind(firstDegree, relative)
      }
      #parents
      else if(!is.null(proband$MotherID) && !is.null(proband$FatherID) && (relative$ID == proband$MotherID || relative$ID == proband$FatherID)){
        firstDegree <-rbind(firstDegree,relative)
      }
    }
  }
  return(firstDegree)
}


#returns all affected family members and proband (affected or otherwise)
affectedFamilyMembers <- function(ped){
  #proband <- ped %>% filter(isProband==1)
  #affectedFam = proband
  for(i in  1:nrow(ped)){
    relative=ped[i,]
    if(relative$ID != proband$ID && relative$isAffAny!=0){
      rbind(affectedFam, proband)
    }

  }
}


first <- firstDegreeFamilyMembers(families[[3]])
View(first)
View(families[[3]])
ambryInfoFam <- ambryRelativeInfo(first)
ambryInfoProband <- ambryProbandInfo(first)
View(ambryInfoFam)
ambryInfo <- rbind.fill(ambryInfoFam, ambryInfoProband)
View(ambryInfo)
affected <- affectedFamilyMembers(ambryInfo)
View(affected)


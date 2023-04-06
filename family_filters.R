#Functions to filter pedigrees in various ways for simulations


removeProbandStatus <- function(ped){
  proband = ped %>% filter(isProband == 1)
  proband$MLH1 = NA
  relatives = ped %>% filter(isProband == 0)
  fam = rbind(proband, relatives)
  return(fam)
}
removeRelativeStatus <- function(ped){
  proband = ped %>% filter(isProband == 1)
  relatives = ped %>% filter(isProband == 0)
  relatives$MLH1 = NA
  fam = rbind(proband, relatives)
  return(fam)
}

firstDegreeFamilyMembers <- function(ped){
  proband <- ped %>% filter(isProband==1)
  firstDegree = proband
  if(!is.null(proband)){
    for(i in  1:nrow(ped)){
      relative=ped[i,]
      #siblings
      if(relative$ID != proband$ID){
        if(!is.null(relative$MotherID) && !is.null(proband$MotherID) && !is.null(relative$FatherID) && !is.null(proband$FatherID) && (relative$MotherID==proband$MotherID) && (relative$FatherID==proband$FatherID)){
          firstDegree <- rbind(firstDegree, relative)
        }
        #parents
        else if(!is.null(proband$MotherID) && !is.null(proband$FatherID) && (relative$ID == proband$MotherID || relative$ID == proband$FatherID)){
          firstDegree <-rbind(firstDegree,relative)
        }
      }
    }
  }
  return(firstDegree)
}


# ambryProbandInfo <- function(ped){
#   proband = ped %>% filter(isProband==1)
#   #proband$riskmod <- character(0)
#   #proband$interAge <- character(0) #should these be character(0)? what does that mean?
#   proband$isDead <- 0
#   return(proband)
# }

### relative info ###
#for simplicity, we will just say that sex can be inferred
#restricting to first degree family
#here genes is just MLH1
# ambryRelativeInfo <- function(ped){
#   relatives = ped %>% filter(isProband==0)
#   #print(relatives)
#   #relatives$CurAge <- NULL
#   relatives$isDead <- NULL
#   #relatives$riskMod <- character(0)
#   #relatives$interAge <- character(0)
#   #View(relatives)
#   return(relatives)
# }


## Function to filter to affected family members, proband, and parents ----
affectedFamilyMembers <- function(ped){
  proband <- ped %>% filter(isProband==1)
  mother <- ped %>% filter(ID == proband$MotherID)
  father <- ped %>% filter(ID == proband$FatherID)
  relatives <- ped %>% filter(isAffAny==1)
  relatives <- relatives %>% filter(isProband==0 && ID != proband$MotherID && ID != proband$FatherID) #to make sure we don't duplicate family members
  affectedFam <- rbind(proband, relatives, mother, father)
  return(affectedFam)
}

### ambry info both proband and relative ---- 
# ambryInfo <- function(ped){
#   proband <- ambryProbandInfo(ped)
#   relatives <- ambryRelativeInfo(ped)
#   family <- rbind.fill(proband, relatives)
#   family <- affectedFamilyMembers(family)
#   return(family)
# }

## Function to mask information from unaffected individuals ----
maskedInfoAmbry <- function(ped){
  #first we will include the proband and the affected relatives
  affected <- ped %>% filter(isProband == 1 | isAffAny == 1)
  #unaffected relatives will have all information masked besides the fact that they exist
  unaffected <- ped %>% filter(isProband == 0 & isAffAny == 0)
  unaffected$CurAge <- 1
  unaffected$isDead <- 0
  unaffected$Twins <- NA
  unaffected$MLH1 <- NA
  masked <- rbind(affected, unaffected)
  return(masked)
}

# maskedMother <- function(ped){
#   #first we will include the proband and the affected relatives
#   pro <- ped %>% filter(isProband == 1)
#   mother <- ped %>% filter(ID == pro$MotherID)
#   other <- ped %>% filter(isProband == 0 && ID != pro$motherID)
#   mother$MLH1 <- NA
#   masked <- rbind(pro, mother)
#   masked <- rbind(masked, other)
#   return(masked)
# }

# affected families lister ----
affectedFamiliesLister <- function(f){
  affectedFamiliesList = list()
  for(i in 1:length(families)){
    fam = as.data.frame(f[[i]])
    fam = fam %>% filter(isAffAny==1) %>% filter(isProband==0)
    if(nrow(fam) > 0){
      affectedFamiliesList[[i]] = nrow(fam)
    }
    else{
      affectedFamiliesList[[i]] = 0}
  }
  return(affectedFamiliesList)
}



# Family Descriptions ----
describeFamilies <- function(fams){
  affectedFamilies = 0 
  affectedProbands = 0
  famSizes = c()
  mlh1Families = 0
  mlh1Probands = 0 
  #affectedFamiliesList = list()
  #j = 1
  for(i in 1:length(fams)){
    f = fams[[i]]
    sizeOfFamily = nrow(f)
    famSizes = c(famSizes, sizeOfFamily)
    f = f %>% filter(isAffAny==1) %>% filter(isProband==0)
    if(nrow(f)>0){
      affectedFamilies = affectedFamilies + 1
      #affectedFamiliesList[[j]] = f
      #j = j + 1
    }
    mlh1fams = fams[[i]] %>% filter(isProband==0) %>% filter(MLH1==1)
    if(nrow(mlh1fams) > 0){
      mlh1Families = mlh1Families + 1
    }
    pb = fams[[i]] %>% filter(isProband==1)
    if(pb$isAffAny == 1){
      affectedProbands = affectedProbands + 1
    }
    if(!is.na(pb$MLH1)){
      if(pb$MLH1 == 1){
        mlh1Probands = mlh1Probands + 1
      }}
  }
  print(paste0("Number of families: ", length(fams)))
  print(paste0("Number of families with affected individuals: ", affectedFamilies))
  print(paste0("Number of families with affected probands: ", affectedProbands))
  print(paste0("Average family size: ", mean(famSizes)))
  print("Summary of family sizes")
  print(summary(famSizes))
  print(paste0("Number of probands that with MLH1 gene: ", mlh1Probands))
  print(paste0("Number of families (non-proband) with MLH1 gene: ", mlh1Families))
}

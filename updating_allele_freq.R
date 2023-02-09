load("~/PanelPRO/data/PanelPRODatabase.rda")
PanelPRODatabase$AlleleFrequency[["MLH1_anyPV",1]] <- 0.5
PanelPRODatabase$AlleleFrequency[["MLH1_anyPV",2]] <- 0.5
PanelPRODatabase$AlleleFrequency[["MLH1_anyPV",3]] <- 0.5
View(PanelPRODatabase$AlleleFrequency)
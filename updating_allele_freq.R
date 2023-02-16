load("~/PanelPRO/data/PanelPRODatabase.rda")
PanelPRODatabase$AlleleFrequency[["MLH1_anyPV",1]] <- 0.2
PanelPRODatabase$AlleleFrequency[["MLH1_anyPV",2]] <- 0.2
PanelPRODatabase$AlleleFrequency[["MLH1_anyPV",3]] <- 0.2
View(PanelPRODatabase$AlleleFrequency)
save(PanelPRODatabase, file = "~/PanelPRO/data/PanelPRODatabase.rda")
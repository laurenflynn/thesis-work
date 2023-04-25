#put panelpro dataframes together
library(tidyverse)
outputs <- "RObjects/summary_tables/job_array_1-400/" #"~/thesis-work/RObjects/panelpro_outputs"

fullTable <- data.frame()
numFiles =0


for(file in list.files(outputs)){
  load(str_interp("RObjects/summary_tables/job_array_1-400/${file}"))
  print(nrow(summaryTable))
  fullTable <- rbind(fullTable, summaryTable)
  numFiles = numFiles  + 1
}

fullTable <- fullTable %>% arrange(famID)
print(nrow(fullTable))
View(fullTable)

summaryTable <- fullTable

save(summaryTable, file = str_interp("RObjects/summary_tables/summaryTable${nrow(fullTable)}Families.Rdata" ))
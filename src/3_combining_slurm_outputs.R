#put panelpro dataframes together
library(tidyverse)
outputs <- "RObjects/summary_tables/no_censoring/"

fullTable <- data.frame()
numFiles =0


for(file in list.files(outputs)){
  load(str_interp("${outputs}/${file}"))
  #print(nrow(summaryTable))
  fullTable <- rbind(fullTable, summaryTable)
  numFiles = numFiles  + 1
}


# create a vector with the expected file numbers
file_nums <- seq(from = 1, to = 10000, by = 25)

# create an empty vector to store the missing file numbers
missing_files <- c()

# loop through the file numbers and check if the file exists
for (i in file_nums) {
  file_name <- paste0("panelPROSummaryTable10000FamiliesNoCensoring", i, "_", i + 24, ".Rdata")
  file_path <- paste0(outputs, file_name)
  if (!file.exists(file_path)) {
    missing_files <- c(missing_files, ((i-1)/25)+1)
  }
}

# print the missing file numbers
cat("Missing files:", paste(missing_files, collapse = ","), "\n")


fullTable <- fullTable %>% arrange(famID)
print(nrow(fullTable))
View(fullTable)

summaryTable <- fullTable

save(summaryTable, file = str_interp("RObjects/summary_tables/summaryTable${nrow(fullTable)}Families_no_censoring.Rdata" ))


#put panelpro dataframes together
library(tidyverse)
outputs <- "RObjects/summary_tables/job_array_1-200/" #"~/thesis-work/RObjects/panelpro_outputs"

fullTable <- data.frame()
numFiles =0


# get a list of all files in the directory
file_list <- list.files(outputs)

# create empty vectors to store the individual numbers
num1 <- numeric()
num2 <- numeric()

# loop through each file name and extract the numbers at the end
for (file_name in file_list) {
  # extract the numbers using a regular expression
  match <- gsub(".*([0-9]+)_([0-9]+)\\.Rdata$", "\\1,\\2", file_name)
  # split the matched string into two separate numbers
  num_vec <- strsplit(match, ",")[[1]]
  num1 <- c(num1, as.numeric(num_vec[1]))
  num2 <- c(num2, as.numeric(num_vec[2]))
}

# print out the unique values in the num1 and num2 vectors
unique_numbers1 <- sort(unique(num1))
unique_numbers2 <- sort(unique(num2))
cat("Numbers found in file names: ", paste(unique_numbers1, unique_numbers2, sep="_", collapse = ", "))

by25 = seq(25,10000,by=25)
by25[!(by25 %in% unique_numbers2)]




for(file in list.files(outputs)){
  load(str_interp("RObjects/summary_tables/job_array_1-200/${file}"))
  print(nrow(summaryTable))
  fullTable <- rbind(fullTable, summaryTable)
  numFiles = numFiles  + 1
}

fullTable <- fullTable %>% arrange(famID)
print(nrow(fullTable))
View(fullTable)

save(fullTable, file = str_interp("~/thesis-work/RObjects/summary_tables/summaryTable${nrow(fullTable)}Families.Rdata" ))
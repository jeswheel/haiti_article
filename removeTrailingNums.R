args <- commandArgs(trailingOnly = TRUE)

removeTrailingNums <- function(file_name) {
  new_name <- gsub("-1-1", '', file_name)
  file.rename(from = file_name, to = new_name)
}

removeTrailingNums(args[1])

args <- commandArgs(trailingOnly = TRUE)

append_bib <- function(file_name) {
  tex_file <- readLines(file_name)

  bib_spot <- which(grepl("\\bibliography{bib-haiti}", tex_file, fixed = TRUE))
  bib_content <- readLines(gsub("\\.tex", "\\.bbl", file_name))

  new_file <- c(tex_file[1:(bib_spot - 1)], bib_content, tex_file[(bib_spot + 1):length(tex_file)])

  writeLines(new_file, con = file_name)
  # file.rename(from = file_name, to = new_name)
}

append_bib(args[1])

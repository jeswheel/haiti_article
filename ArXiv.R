args <- commandArgs(trailingOnly = TRUE)

ms_to_ArXiv <- function(file_name) {
  tmp_file <- readLines(file_name)
  tmp_file <- gsub("figure/", "", tmp_file)
  tmp_file <- gsub("\\linenumbers", "", tmp_file, fixed = TRUE)
  tmp_file <- gsub("\\nolinenumber", "", tmp_file, fixed = TRUE)
  lineno_pkg <- which(grepl('lineno', tmp_file))
  tmp_file <- c(tmp_file[1:(lineno_pkg - 1)], tmp_file[(lineno_pkg + 1):length(tmp_file)])

  # geometry_line <- which(grepl("{geometry}", tmp_file, fixed = TRUE))
  # tmp_file[geometry_line] <- "\\usepackage[top=0.85in,left=1.5in,footskip=0.75in]{geometry}"
  # tmp_file <- gsub("bib-haiti", 'ms', tmp_file)

  bib_spot <- which(grepl("\\bibliography{bib-haiti,bib-lit}", tmp_file, fixed = TRUE))
  bib_content <- readLines(gsub("\\.tex", "\\.bbl", file_name))

  new_file <- c(tmp_file[1:(bib_spot - 1)], bib_content, tmp_file[(bib_spot + 1):length(tmp_file)])

  writeLines(new_file, con = file_name)
}

si_to_ArXiv <- function() {
  tmp_file <- readLines("ArXiv/si.tex")
  tmp_file <- gsub("inputs/", "", tmp_file)
  tmp_file <- gsub("figure/", "", tmp_file)
  # tmp_file <- gsub("../bib-haiti", 'ms', tmp_file)

  input_lines <- which(grepl("\\\\input\\{([^}]*)\\}", tmp_file))
  input_files <- paste0(gsub("\\\\input\\{|\\}", "", tmp_file[input_lines]), ".tex")

  for (i in 1:length(input_files)) {
    tmp_input_file <- readLines(paste0("ArXiv/", input_files[i]))
    tmp_input_file <- gsub("figure/", "", tmp_input_file)
    writeLines(tmp_input_file, con = paste0("ArXiv/", input_files[i]))
  }

  bib_spot <- which(grepl("bib-haiti", tmp_file, fixed = TRUE))
  bib_content <- readLines("ArXiv/si.bbl")

  new_file <- c(tmp_file[1:(bib_spot - 1)], bib_content, tmp_file[(bib_spot + 1):length(tmp_file)])

  writeLines(new_file, con = "ArXiv/si.tex")
}

ms_to_ArXiv(args[1])
si_to_ArXiv()

remove_start_end <- function(file) {
  tmp_file <- readLines(file)
  start_line <- which(grepl("%* START", tmp_file))
  end_line <- which(grepl("%* END", tmp_file))
  tmp_file[start_line:end_line]
  writeLines(tmp_file[(start_line+1):(end_line-1)], con = gsub("\\.tex", "Out.tex", file))
}

make_si_file <- function(file) {
  template <- readLines('inputs/si_files_template.tex')
  input_file <- readLines(file)

  if (grepl('Tab', file)) {
    tab_caption <- which(grepl('\\caption', input_file))
    input_file <- input_file[-tab_caption]
  }

  start_line <- which(grepl("%* START", template))
  end_line <- which(grepl("%* END", template))
  file_loc <- gsub("inputs/", "", file)
  new_file <- c(template[1:(start_line - 1)], input_file, template[(end_line + 1):length(template)])
  writeLines(new_file, con = file_loc)

  has_citation <- any(grepl("\\cite", new_file))
  if (has_citation) writeLines(as.character(has_citation), con = gsub('.tex', '.tmp', file_loc))
}

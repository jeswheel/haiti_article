remove_start_end <- function(file) {
  tmp_file <- readLines(file)
  start_line <- which(grepl("%* START", tmp_file))
  end_line <- which(grepl("%* END", tmp_file))
  tmp_file[start_line:end_line]
  writeLines(tmp_file[(start_line+1):(end_line-1)], con = gsub("\\.tex", "Out.tex", file))
}


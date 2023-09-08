args <- commandArgs(trailingOnly = TRUE)

remove_knitrout <- function(file_name) {
  tmp_file <- readLines(file_name)
  knitrout_env <- which(grepl("\\newenvironment{knitrout}{}{}", tmp_file, fixed = TRUE))
  begin_env <- which(grepl("\\begin{knitrout}", tmp_file, fixed = TRUE))
  end_env <- which(grepl("\\end{knitrout}", tmp_file, fixed = TRUE))

  if (length(begin_env) != length(end_env)) {
    stop("Number of \\begin{knitrout} does not match number of \\end{knitrout}.")
  }

  row_ids <- c(knitrout_env)

  for (i in 1:length(begin_env)) {
    if (begin_env[i] > end_env[i]) {
      stop("Problem matching \\begin{knitrout} with \\end{knitrout}")
    }

    row_ids <- c(row_ids, begin_env[i]:end_env[i])
  }

  writeLines(tmp_file[-row_ids], con = gsub("\\.tex", "-submission.tex", file_name))
  # file.rename(from = file_name, to = new_name)
}

remove_knitrout(args[1])

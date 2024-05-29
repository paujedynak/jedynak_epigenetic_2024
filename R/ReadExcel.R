#' Read several excel tabs at once
#'
#' @param filename A string defining an excel file to be read
#' @param tibble Default FALSE
#'
#' @return A list where each sheet is a list element
#' @export
#' @import readxl

read_excel_allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- excel_sheets(filename)
  x <- lapply(sheets, function(X) read_excel(filename, sheet = X))
  if (!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}
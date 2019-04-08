#' User inputs
#'
#' Function with no parameters, requests inputs from user to run enrichment.
user_inputs <- function() {
  selection <- readline(promp = "Number of random sets: ")
  selection <- as.integer(selection)
  return(selection)
}

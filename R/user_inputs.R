#' User inputs
#'
#' Function with no parameters, requests inputs from user to run enrichment.
user_inputs <- function() {

  fthr <- readline(promp = "Fold change threshold: ")
  if (fthr == "") {
    message("No fold change threshold selected. Setting to default (F = 0).")
    fthr <- 0
  } else {
    fthr <- abs(as.numeric(fthr))
  }

  pthr <- readline(promp = "p-value threshold: ")
  if (pthr == "") {
    message("No p-value threshold selected. Setting to default (p < 0.05).")
    pthr <- 0.05
  } else if (pthr > 1) {
    message("Invalid p-value threshold. Setting to default (p < 0.05).")
    pthr <- 0.05
  } else {
    pthr <- as.numeric(pthr)
  }

  nr <- readline(promp = "Number of random sets: ")
  if (nr == "") {
    message("No number of random sets selected. Setting to default (N = 1000).")
    nr <- 1000
  } else {
    nr <- as.integer(nr)
  }


  ui <- list(fthr, pthr, nr)
  return(ui)
}

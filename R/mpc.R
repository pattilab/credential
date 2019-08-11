#' Minimum and maximum mass-per-carbon values of metabolites in E.coli 
#'
#' The minimum and maximum mass-per-carbons are used to refine the isotopologue searching.  
#' 
#' @format A data.frame with three variables and 320 entries.
#' \describe{
#' \item{carbons}{The number of carbons}
#' \item{min_mpc}{maximum mass-per-carbon among metabolites with designated carbon number}
#' \item{max_mpc}{minimum mass-per-carbon among metabolites with designated carbon number}
#' }

#' @source \url{http://ecmdb.ca/}
"mpc"
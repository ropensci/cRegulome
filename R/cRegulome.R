#' \code{cRegulome} package
#'
#' Download, access and visualize Regulome (microRNA and transcription factors)
#' data from miRCancer and Cistrome cancer
#'
#' @docType package
#' @name cRegulome
#' 
#' @import DBI
#' @import RSQLite
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
## fix by @jennybc
## source https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
if(getRversion() >= "2.15.1")  utils::globalVariables(c('study',
                                                        'mirna_base',
                                                        'cor',
                                                        'tf',
                                                        'feature'))

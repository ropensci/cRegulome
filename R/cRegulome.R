#' \code{cRegulome} package
#'
#' Download, access and visualize Regulome (microRNA and transcription factors)
#' data from miRCancer and Cistrome cancer
#' 
#' @section \code{cRegulome} functions to download and query the database file:
#' \code{\link{get_db}}
#' \code{\link{get_tf}}
#' \code{\link{get_mir}}
#' 
#' @section \code{cRegulome} functions to create S3 objects:
#' \code{\link{cTF}}
#' \code{\link{cmicroRNA}}
#' 
#' @section \code{cRegulome} functions to reshape S3 objects:
#' \code{\link{cor_tidy}}
#' \code{\link{cor_igraph}}
#' 
#' @section \code{cRegulome} functions to visualize data in S3 objects:
#' \code{\link{cor_hist}}
#' \code{\link{cor_joy}}
#' \code{\link{cor_plot}}
#' \code{\link{cor_upset}}
#' \code{\link{cor_venn_diagram}}
#' 
#' @docType package
#' @name cRegulome
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
## fix by @jennybc
## source https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
#if(getRversion() >= "2.15.1")  utils::globalVariables(c('study',
#                                                        'mirna_base',
#                                                        'cor',
#                                                        'tf',
#                                                        'feature'))

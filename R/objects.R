#' Construct cmicroRNA object
#'
#' Constructs an S3 object called cmicroRNA contains data returned by calling
#' \link{get_mir}. Used to define methods for printing and visualizing
#' microRNA-gene expression correlations.
#'
#' @param dat_mir A \code{data.frame} such as this returned by calling
#'    \link{get_mir}.
#'
#' @return An S3 object of class \code{cmicroRNA}
#' @examples
#' # connect to test database file
#' db_file <- system.file("extdata", "cRegulome.db", package = "cRegulome")
#' conn <- DBI::dbConnect(RSQLite::SQLite(), db_file)
#'
#' # optional arguments
#' dat <- get_mir(conn,
#'     mir = c('hsa-let-7b', 'hsa-mir-134'),
#'     study = c('ACC', 'BLCA'),
#'     min_cor = .5,
#'     max_num = 100,
#'     targets_only = TRUE)
#' DBI::dbDisconnect(conn)
#'
#' # convert object
#' ob <- cmicroRNA(dat)
#'
#' @export
cmicroRNA <- function(dat_mir){
    microRNA <- unique(dat_mir$mirna_base)
    features <- unique(dat_mir$feature)
    studies <- unique(dat_mir$study)

    `%>%` <- dplyr::`%>%`

    if(length(studies) == 1) {
        corr <- dat_mir %>%
            reshape2::dcast(feature ~ mirna_base, value.var = 'cor')
    } else {
        corr <- lapply(unique(dat_mir$study),
                       function(x) {
                           dat_mir %>%
                               dplyr::filter(study == x) %>%
                               reshape2::dcast(feature ~ mirna_base,
                                               value.var = 'cor')
                       })
        names(corr) <- unique(dat_mir$study)
    }
    structure(list(
        microRNA = microRNA,
        features = features,
        studies = studies,
        corr = corr
    ),
    class = 'cmicroRNA')
}

#' Construct cTF object
#'
#' Constructs an S3 object called cTF contains data returned by calling
#' \link{get_tf}. Used to define methods for printing and visualizing
#' transcription factors-gene expression correlations.
#'
#' @param dat_tf A \code{data.frame} such as this returned by calling
#'    \link{get_tf}.
#'
#' @return An S3 object of class \code{cTF}
#' @examples
#' # connect to test database file
#' db_file <- system.file("extdata", "cRegulome.db", package = "cRegulome")
#' conn <- DBI::dbConnect(RSQLite::SQLite(), db_file)
#'
#' # optional arguments
#' dat <- get_tf(conn,
#'     tf = c('AFF4', 'ESR1'),
#'     study = c('ACC', 'BLCA'),
#'     min_cor = .5,
#'     max_num = 100,
#'     targets_only = TRUE)
#' DBI::dbDisconnect(conn)
#'
#' # convert to object
#' ob <- cTF(dat)
#'
#' @export
cTF <- function(dat_tf){
    TF <- unique(dat_tf$tf)
    features <- unique(dat_tf$feature)
    studies <- unique(dat_tf$study)

    `%>%` <- dplyr::`%>%`

    if(length(studies) == 1) {
        corr <- dat_tf %>%
            reshape2::dcast(feature ~ tf, value.var = 'cor')
    } else {
        corr <- lapply(unique(dat_tf$study),
                       function(x) {
                           dat_tf %>%
                               dplyr::filter(study == x) %>%
                               reshape2::dcast(feature ~ tf, value.var = 'cor')
                       })
        names(corr) <- unique(dat_tf$study)
    }
    structure(list(
        TF = TF,
        features = features,
        studies = studies,
        corr = corr
    ),
    class = 'cTF')
}

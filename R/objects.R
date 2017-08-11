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
#' # downlaod db file
#' get_db(test = TRUE)
#'
#' # get data for 2 microRNAs in the ACC study
#' dat <- get_mir(c('hsa-let-7b', 'hsa-mir-134'),
#'     study = c('ACC', 'BLCA'),
#'     min_cor = .5)
#'
#' # convert to cmicroRNA object
#' ob <- cmicroRNA(dat)
#'
#' @import dplyr reshape2
#' @export
cmicroRNA <- function(dat_mir){
    microRNA <- unique(dat_mir$mirna_base)
    features <- unique(dat_mir$feature)
    studies <- unique(dat_mir$study)
    if(length(studies) == 1) {
        corr <- dat_mir %>%
            dcast(feature ~ mirna_base, value.var = 'cor')
    } else {
        corr <- lapply(unique(dat_mir$study),
                       function(x) {
                           dat_mir %>%
                               filter(study == x) %>%
                               dcast(feature ~ mirna_base, value.var = 'cor')
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
#' # downlaod db file
#' get_db(test = TRUE)
#'
#' # get data for 2 tarnscription factors in the ACC study
#' dat <- get_tf(c('AFF4', 'ESR1'),
#'     study = c('ACC', 'BLCA'),
#'     min_cor = .5)
#'
#' # convert to cTF object
#' ob <- cTF(dat)
#'
#' @import dplyr reshape2
#' @export
cTF <- function(dat_tf){
    TF <- unique(dat_tf$tf)
    features <- unique(dat_tf$feature)
    studies <- unique(dat_tf$study)
    if(length(studies) == 1) {
        corr <- dat_tf %>%
            dcast(feature ~ tf, value.var = 'cor')
    } else {
        corr <- lapply(unique(dat_tf$study),
                       function(x) {
                           dat_tf %>%
                               filter(study == x) %>%
                               dcast(feature ~ tf, value.var = 'cor')
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

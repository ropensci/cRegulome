#' Construct cmicroRNA object
#'
#' Constructs an S3 object called cmicroRNA contains data returned by calling
#' \link{get_mir}. Used to define methods for printing and visualizing
#' microRNA-gene expression correlations.
#'
#' @param dat_mir A \code{data.frame} such as this returned by calling
#' \link{get_mir}.
#'
#' @return An S3 object of class \code{cmicroRNA}
#' 
#' @examples 
#' # load required libraries
#' library(RSQLite)
#' library(cRegulome)
#' 
#' # locate the testset file and connect
#' fl <- system.file('extdata', 'cRegulome.db', package = 'cRegulome')
#' conn <- dbConnect(SQLite(), fl)
#' 
#' # enter a custom query with different arguments
#' dat <- get_mir(conn,
#'                mir = 'hsa-let-7g',
#'                study = 'STES',
#'                min_abs_cor = .3,
#'                max_num = 5)
#' 
#' # make a cmicroRNA object   
#' cmir <- cmicroRNA(dat)
#' 
#' @importFrom magrittr %>%
#' 
#' @export
cmicroRNA <- function(dat_mir){
    # extract items of the list from the data.frame
    microRNA <- unique(dat_mir$mirna_base)
    features <- unique(dat_mir$feature)
    studies <- unique(dat_mir$study)


        
    # reshape the data.frame/s 
    # microRNA in columns and feature in rows
    # object contains a single study
    if(length(studies) == 1) {
        corr <- dat_mir %>%
            reshape2::dcast(feature ~ mirna_base, value.var = 'cor')
    } else {
        # object with multiple studies
        # returns a list of data.frames
        corr <- lapply(unique(dat_mir$study),
                    function(x) {
                        dat_mir %>%
                            dplyr::filter(study == x) %>%
                            reshape2::dcast(feature ~ mirna_base,
                                            value.var = 'cor')
                        })
        names(corr) <- unique(dat_mir$study)
    }
    
    # object structure and class name
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
#' \link{get_tf}.
#'
#' @return An S3 object of class \code{cTF}
#' 
#' @examples 
#' # load required libraries
#' library(RSQLite)
#' library(cRegulome)
#' 
#' # locate the testset file and connect
#' fl <- system.file('extdata', 'cRegulome.db', package = 'cRegulome')
#' conn <- dbConnect(SQLite(), fl)
#' 
#' # enter a custom query with different arguments
#' dat <- get_tf(conn,
#'               tf = 'LEF1',
#'               study = 'STES*',
#'               min_abs_cor = .3,
#'               max_num = 5)
#' 
#' # make a cTF object   
#' ctf <- cTF(dat)
#' 
#' @importFrom magrittr %>%
#' 
#' @export
cTF <- function(dat_tf){
    # extract items of the list from the data.frame
    TF <- unique(dat_tf$tf)
    features <- unique(dat_tf$feature)
    studies <- unique(dat_tf$study)


        
    # reshape the data.frame/s 
    # tf in columns and feature in rows
    # object contains a single study
    if(length(studies) == 1) {
        corr <- dat_tf %>%
            reshape2::dcast(feature ~ tf, value.var = 'cor')
    } else {
        # object with multiple studies
        # returns a list of data.frames
        corr <- lapply(unique(dat_tf$study),
                    function(x) {
                        dat_tf %>%
                            dplyr::filter(study == x) %>%
                            reshape2::dcast(feature ~ tf, value.var = 'cor')
                        })
        names(corr) <- unique(dat_tf$study)
    }
    # object structure and class name
    structure(list(
        TF = TF,
        features = features,
        studies = studies,
        corr = corr
    ),
    class = 'cTF')
}

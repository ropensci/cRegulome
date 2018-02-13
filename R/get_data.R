#' Get cRegulome.db file
#'
#' This function calls \code{\link[utils]{download.file}} to download the
#' pre-build database file of cRegulome. Additionally, the function checks
#' the validity of the pre-defined URL and whether the database file exists
#' in the current working directory to avoid redownloading it. Typically,
#' users would run this function once at the first time the use the package
#' or to update the database to the latest version.
#'
#' @param test A \code{logical}, default \code{FALSE}. When \code{TRUE}
#' downloads a database file with the same structure with a subset of
#' the data for speed.
#' @param destfile A character vector for the desired path for the database 
#' file. By default, when not specified, is constructed by using 
#' \code{\link{tempdir}} as a directory and the string \code{cRegulome.db.gz}
#' @param ... Optional arguments passed to \code{\link[utils]{download.file}}
#'
#' @return Downloads a compressed \code{sqlite} file to the current working
#' directory. The file is named \code{cRegulome.db.gz} by default and it's
#' not advised to change the name to avoid breaking the other functions
#' that calls the database.
#'
#' @examples 
#' \dontrun{
#' # download a test set of the database
#' get_db(test = TRUE)
#' 
#' # download the full database file
#' get_db(test = FALSE)
#' }
#' 
#'
#' # load the test db file from shipped with the pacakge
#' db_file <- system.file("extdata", "cRegulome.db", package = "cRegulome")
#' file.info(db_file)
#'
#' @export
get_db <- function(test = FALSE, destfile, ...) {
    # db file url
    if(test == TRUE) {
        url <- 'https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/10330329/test.db.gz'
    } else {
        url <- 'https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/9537385/cRegulome.db.gz'
    }

    # check url exists
    if(httr::http_error(url)) {
        stop("URL doesn't exist.")
    }
    
    # make a destfile
    if(missing(destfile)) {
        destfile <- paste(tempdir(), 'cRegulome.db.gz', sep = '/')
    }
    # download file
    if(file.exists(destfile)) {
        message('File already exists in the directory.')
    } else {
        tryCatch({
            utils::download.file(url, 
                                 destfile = destfile,
                                 mode = 'wb')
            R.utils::gunzip(destfile)
            },
                error = function(){
                    message('Downloading file failed')
                    return(NA)
                    })
    }
    }

#' Get microRNA correlations from cRegulome.db
#'
#' This function access the \code{sqlite} database file which is obtained by
#' running \link{get_db}. Basically, the function provides ways to query the 
#' database to the correlation data of the microRNAs of interest. The function 
#' returns an error if the database file \code{cRegulome.db} is not in the 
#' working directory.
#'
#' @param conn A connection to the database file by \code{\link[DBI]{dbConnect}}
#' @param mir A required \code{character} vector of the microRNAs of interest.
#' These are the miRBase ID which are the official identifiers of the
#' widely used miRBase database, \url{http://www.mirbase.org/}.
#' @param study A \code{character} vector of The Cancer Genome Atlas (TCGA)
#' study identifiers. To view the available studies in TCGA project,
#' \url{https://tcga-data.nci.nih.gov/docs/publications/tcga}. When left to
#' default \code{NULL} all available studies will be included.
#' @param min_abs_cor A \code{numeric}, an absolute correlation minimum between 0
#' and 1 for each \code{mir}.
#' @param max_num An \code{integer}, maximum number of \code{features} to show
#' for each \code{mir} in each \code{study}.
#' @param targets_only A \code{logical}, default \code{FALSE}. When
#' \code{TRUE}, \code{features} will be the microRNA targets as defined in
#' the package targetscan.Hs.eg.db.
#'
#' @return A tidy \code{data.frame} of four columns. \code{mirna_base} is the
#' microRNA miRBase IDs, \code{feature} is the features/genes, \code{cor}
#' is the corresponding expression correlations and \code{study} is TCGA
#' study ID.
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
#' # get microRNA correlations in all studies
#' get_mir(conn,
#'         mir = 'hsa-let-7g')
#' 
#' # get correlations in a particular study
#' get_mir(conn,
#'         mir = 'hsa-let-7g',
#'         study = 'STES')
#' 
#' # enter a custom query with different arguments
#' get_mir(conn,
#'         mir = 'hsa-let-7g',
#'         study = 'STES',
#'         min_abs_cor = .3,
#'         max_num = 5)
#'         
#' @importFrom magrittr %>%
#' 
#' @export
get_mir <- function(conn,
                    mir,
                    study = NULL,
                    min_abs_cor = NULL,
                    max_num = NULL,
                    targets_only = FALSE) {

    # unpack filters and check types
    table <- 'cor_mir'

    if(is.null(mir)) {
        stop("User should provide at least one microRNA ID")
    } else if(!is.character(mir)) {
        stop("mir should be a character vector")
    } else {
        mir <- as.character(mir)
    }

    if(is.null(study)) {
        study <- DBI::dbListFields(conn, table)[-1:-2]
    } else if(!is.character(study)){
        stop("Study should be a character vector")
    } else {
        study <- as.character(study)
    }

    if(is.null(min_abs_cor)) {
        min_abs_cor <- 0
    } else if(!is.numeric(min_abs_cor) || min_abs_cor > 1 || min_abs_cor < 0) {
        stop("min_abs_cor should be a numeric between 0 and 1.")
    }
    else {
        min_abs_cor <- as.numeric(min_abs_cor)
    }

    # get main data by applying filters and tidy
    
    dat <- conn %>%
        dplyr::tbl(table) %>%
        dplyr::select(mirna_base, feature, study) %>%
        dplyr::filter(mirna_base %in% mir) %>%
        dplyr::collect() %>%
        tidyr::gather(study, cor, -mirna_base, -feature) %>%
        dplyr::mutate(cor = cor/100) %>%
        dplyr::filter(abs(cor) > min_abs_cor) %>%
        dplyr::arrange(dplyr::desc(abs(cor))) %>%
        stats::na.omit()

    # apply targets only filters when TRUE
    if(targets_only == TRUE) {
        # subset targets
        targets <- conn %>%
            dplyr::tbl('targets_mir') %>%
            dplyr::filter(mirna_base %in% mir) %>%
            dplyr::collect() %>%
            unique()

        # subset main data to targets only
        dat <- dplyr::inner_join(dat, targets) %>%
            stats::na.omit()
    }

    # subset to max_num
    if(is.null(max_num)) {
    } else if(max_num %% 1 != 0 || max_num <= 0) {
        stop("max_num should be a positive integer.")
    } else {
        dat <- dat %>%
            dplyr::group_by(mirna_base, study) %>%
            dplyr::slice(1:max_num)
    }

    # return dat
    return(dat)
}

#' Get transcription factor correlations from cRegulome.db
#'
#' This function access the \code{sqlite} database file which is obtained by
#' running \link{get_db}. Basically, the function provides ways to query the 
#' database to the correlation data of the transcription factors of interest. 
#' The function returns an error if the database file \code{cRegulome.db} is 
#' not in the working directory.
#'
#' @param tf A required \code{character} vector of the transcription factor of
#' interest. These are the HUGO official gene symbols of the genes contains the
#' transcription factor.
#' @inheritParams get_mir
#' @param min_abs_cor A \code{numeric}, an absolute correlation minimum between
#'  0 and 1 for each \code{tf}.
#' @param max_num An \code{integer}, maximum number of \code{features} to show
#' for each \code{tf} in each \code{study}.
#' @param targets_only A \code{logical}, default \code{FALSE}. When
#' \code{TRUE}, \code{features} will be the targets of the transcription
#' factors as identified in the Cistrome Cancer,
#' \url{http://cistrome.org/CistromeCancer/}
#'
#' @return A tidy \code{data.frame} of four columns. \code{tf} is the official
#' gene symbols of the genes contains the transcription factor, \code{feature}
#' is the features/genes, cor is the corresponding expression correlations
#' and \code{study} is TCGA study ID.
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
#' # get transcription factors correlations in all studies
#' get_tf(conn,
#'         tf = 'LEF1')
#' 
#' # get correlations in a particular study
#' get_tf(conn,
#'        tf = 'LEF1',
#'        study = 'STES*')
#' 
#' # enter a custom query with different arguments
#' get_tf(conn,
#'        tf = 'LEF1',
#'        study = 'STES*',
#'        min_abs_cor = .3,
#'        max_num = 5)
#'        
#' @importFrom magrittr %>%
#' 
#' @export
get_tf <- function(conn,
                    tf,
                    study = NULL,
                    min_abs_cor = NULL,
                    max_num = NULL,
                    targets_only = FALSE) {

    # unpack filters and check types
    table <- 'cor_tf'

    if(is.null(tf)) {
        stop("User should provide at least one TF ID")
    } else if(!is.character(tf)) {
        stop("tf should be a character vector")
    } else {
        tf_id <- as.character(tf)
    }

    if(is.null(study)) {
        studies <- DBI::dbListFields(conn, table)[-1:-2]
    } else if(!is.character(study)){
        stop("Study should be a character vector")
    } else {
        studies <- as.character(study)
    }

    if(is.null(min_abs_cor)) {
        min_abs_cor <- 0
    } else if(!is.numeric(min_abs_cor) || min_abs_cor > 1 || min_abs_cor < 0) {
        stop("min_abs_cor should be a numeric between 0 and 1.")
    }
    else {
        min_abs_cor <- as.numeric(min_abs_cor)
    }

    # get main data by applying filters and tidy
    
    dat <- conn %>%
        dplyr::tbl(table) %>%
        dplyr::select(tf, feature, studies) %>%
        dplyr::filter(tf %in% tf_id) %>%
        dplyr::collect() %>%
        tidyr::gather(study, cor, -tf, -feature) %>%
        dplyr::mutate(cor = cor/100) %>%
        dplyr::filter(abs(cor) > min_abs_cor) %>%
        dplyr::arrange(dplyr::desc(abs(cor))) %>%
        stats::na.omit()

    # apply targets only filters when TRUE
    if(targets_only == TRUE) {
        # subset targets
        tf_id <- unique(dat$tf)
        targets <- conn %>%
            dplyr::tbl('targets_tf') %>%
            dplyr::filter(tf %in% tf_id, study %in% studies) %>%
            dplyr::collect()

        # subset main data to targets only
        dat <- dplyr::inner_join(dat, targets) %>%
            stats::na.omit()
    }

    # subset to max_num
    if(is.null(max_num)) {
    } else if(max_num %% 1 != 0 | max_num <= 0) {
        stop("max_num should be a positive integer.")
    } else {
        dat <- dat %>%
            dplyr::group_by(tf, study) %>%
            dplyr::slice(1:max_num)
    }

    # return dat
    return(dat)
}

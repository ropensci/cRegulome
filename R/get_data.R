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
#'    downlaods a database file with the same structure with a subset of
#'    the data for speed.
#' @param ... Optional arguments passed to \code{\link[utils]{download.file}}
#'
#' @return Downloads a compressed \code{sqlite} file to the current working
#'    directory. The file is named \code{cRegulome.db.gz} by default and it's
#'    not advised to change the name to avoid breaking the other functions
#'    that calls the database.
#' @examples
#' \dontrun{
#' # downlaod db file
#' get_db(test = TRUE)
#'
#' # check it exits in the current working directory
#' # should return TRUE
#' file.exists('./cRegulome.db.gz')
#' }
#' @export
get_db <- function(test = FALSE, ...) {
    # db file url
    if(test == TRUE) {
        url <- 'https://www.dropbox.com/s/t8ga5j8o81jkcuv/test.db.gz?raw=1'
    } else {
        url <- 'https://www.dropbox.com/s/t8ga5j8o81jkcuv/cRegulome.db.gz?raw=1'
    }

    # check url exists
    if(!RCurl::url.exists(url)) {
        stop("URL doesn't exist.")
    }

    # download file
    if(file.exists('cRegulome.db')) {
        message('File already exists in the current directory.')
    } else {
        tryCatch(utils::download.file(url, destfile = 'cRegulome.db.gz'),
                 error = function(){
                     message('Downloading file failed')
                     return(NA)
                 })
    }
    }

#' Get microRNA correlations from cRegulome.db
#'
#' This function access the \code{sqlite} database file which is obtained by
#' running \link{get_db}.The function returns an error if the uncompressd
#' database file \code{cRegulome.db} is not in the working directory.
#' Basically, the function provides fileters to subset the large tables to
#' the items of interest.
#'
#' @param conn A connection to the database file by \code{\link[DBI]{dbConnect}}
#' @param mir A required \code{character} vector of the microRNAs of interest.
#'    These are the miRBase ID which are the official identifiers of the
#'    widely used miRBase database, \url{http://www.mirbase.org/}.
#' @param study A \code{character} vector of The Cancer Genome Atlase (TCGA)
#'    study identifiers. To view the available studies in TCGA project,
#'    \url{https://tcga-data.nci.nih.gov/docs/publications/tcga}. When left to
#'     defult \code{NULL} all available studies will be included.
#' @param min_cor A \code{numeric}, an absoute correlation minimmum between 0
#'     and 1 for each \code{mir}.
#' @param max_num An \code{integer}, maximum number of \code{features} to show
#'     for each \code{mir} in each \code{study}.
#' @param targets_only A \code{logical}, default \code{FALSE}. When
#'    \code{TRUE}, \code{features} will be the microRNA targets as defined in
#'    the package \code{\link[targetscan.Hs.eg.db]{targetscan.Hs.eg.db}}.
#'
#' @return A tidy \code{data.frame} of four columns. \code{mirna_base} is the
#'    microRNA miRBase IDs, \code{feature} is the features/genes, \code{cor}
#'    is the corresponding expression correaltions and \code{study} is TCGA
#'    study ID.
#' @examples
#' \dontrun{
#' # downlaod db file
#' get_db(test = TRUE)
#' gunzip('cRegulome.db.gz')
#' }
#' # connect to test database file
#' db_file <- system.file("extdata", "cRegulome.db", package = "cRegulome")
#' conn <- DBI::dbConnect(RSQLite::SQLite(), db_file)
#'
#' # provide only required arguments
#' dat <- get_mir(conn, mir = 'hsa-let-7b')
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
#' @export
get_mir <- function(conn,
                    mir,
                    study = NULL,
                    min_cor = NULL,
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

    if(is.null(min_cor)) {
        min_cor <- 0
    } else if(!is.numeric(min_cor) || min_cor > 1 || min_cor < 0) {
        stop("min_cor should be a numeric between 0 and 1.")
    }
    else {
        min_cor <- as.numeric(min_cor)
    }

    # get main data by applying filters and tidy
    `%>%` <- dplyr::`%>%`

    dat <- conn %>%
        dplyr::tbl(table) %>%
        dplyr::select(mirna_base, feature, study) %>%
        dplyr::filter(mirna_base %in% mir) %>%
        dplyr::collect() %>%
        tidyr::gather(study, cor, -mirna_base, -feature) %>%
        dplyr::mutate(cor = cor/100) %>%
        dplyr::filter(abs(cor) > min_cor) %>%
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
    } else if(is.integer(max_num) || max_num < 0) {
        stop("max_num should be an integer between 1 and Inf.")
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
#' running \link{get_db}. The function returns an error if the uncompressd
#' database file \code{cRegulome.db} is not in the working directory.
#' Basically, the function provides fileters to subset the large tables to the
#' items of interest.
#'
#' @param tf A required \code{character} vector of the transcription factor of
#'    interest. These are the official gene symbols of the genes contains the
#'    transcription factor.
#' @inheritParams get_mir
#' @param min_cor A \code{numeric}, an absoute correlation minimmum between 0
#'    and 1 for each \code{tf}.
#' @param max_num An \code{integer}, maximum number of \code{features} to show
#'    for each \code{tf} in each \code{study}.
#' @param targets_only A \code{logical}, default \code{FALSE}. When
#'    \code{TRUE}, \code{features} will be the targets of the transcription
#'    factors as identified in the Cistrome Cancer,
#'    \url{http://cistrome.org/CistromeCancer/}
#'
#' @return A tidy \code{data.frame} of four columns. \code{tf} is the official
#'    gene symbols of the genes contains the transcription factor,
#'    \code{feature} is the features/genes, cor is the corresponding
#'    expression correaltions and \code{study} is TCGA study ID.
#' @examples
#' \dontrun{
#' # downlaod db file
#' get_db(test = TRUE)
#' R.utils::gunzip('cRegulome.db.gz')
#' }
#' # connect to test database file
#' db_file <- system.file("extdata", "cRegulome.db", package = "cRegulome")
#' conn <- DBI::dbConnect(RSQLite::SQLite(), db_file)
#'
#' # provide only required arguments
#' dat <- get_tf(conn, tf = 'AFF4')
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
#' @export
get_tf <- function(conn,
                   tf,
                   study = NULL,
                   min_cor = NULL,
                   max_num = NULL,
                   targets_only = FALSE) {

    # unpack filters and check types
    table <- 'cor_tf'

    if(is.null(tf)) {
        stop("User should provide at least one microRNA ID")
    } else if(!is.character(tf)) {
        stop("tf should be a character vector")
    } else {
        tf_id <- as.character(tf)
    }

    if(is.null(study)) {
        study <- DBI::dbListFields(conn, table)[-1:-2]
    } else if(!is.character(study)){
        stop("Study should be a character vector")
    } else {
        study <- as.character(study)
    }

    if(is.null(min_cor)) {
        min_cor <- 0
    } else if(!is.numeric(min_cor) || min_cor > 1 || min_cor < 0) {
        stop("min_cor should be a numeric between 0 and 1.")
    }
    else {
        min_cor <- as.numeric(min_cor)
    }

    # get main data by applying filters and tidy
    `%>%` <- dplyr::`%>%`
    dat <- conn %>%
        dplyr::tbl(table) %>%
        dplyr::select(tf, feature, study) %>%
        dplyr::filter(tf %in% tf_id) %>%
        dplyr::collect() %>%
        tidyr::gather(study, cor, -tf, -feature) %>%
        dplyr::mutate(cor = cor/100) %>%
        dplyr::filter(abs(cor) > min_cor) %>%
        dplyr::arrange(desc(abs(cor))) %>%
        stats::na.omit()

    # apply targets only filters when TRUE
    if(targets_only == TRUE) {
        # subset targets
        tf_id <- unique(dat$tf)
        targets <- conn %>%
            dplyr::tbl('targets_tf') %>%
            dplyr::filter(tf %in% tf_id) %>%
            dplyr::collect()

        # subset main data to targets only
        dat <- dplyr::inner_join(dat, targets) %>%
            stats::na.omit()
    }

    # subset to max_num
    if(is.null(max_num)) {
    } else if(is.integer(max_num) || max_num < 0) {
        stop("max_num should be an integer between 1 and Inf.")
    } else {
        dat <- dat %>%
            dplyr::group_by(tf, study) %>%
            dplyr::slice(1:max_num)
    }

    # return dat
    return(dat)
}

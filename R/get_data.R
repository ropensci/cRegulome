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
#' # load the test db file from shipped with the pacakge
#' db_file <- system.file("extdata", "cRegulome.db", package = "cRegulome")
#' file.info(db_file)
#' 
#' @importFrom httr http_error
#' @importFrom utils download.file
#' @importFrom R.utils gunzip
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
    if(http_error(url)) {
        stop("URL doesn't exist.")
    }
    
    # make a destfile
    if(missing(destfile)) {
        destfile <- paste(tempdir(), 'cRegulome.db', sep = '/')
    }
    # download file
    if(file.exists(destfile)) {
        message('File already exists in the directory.')
    } else {
        download.file(url, 
                      destfile = paste(destfile, 'gz', sep = '.'),
                      mode = 'wb')
        gunzip(paste(destfile, 'gz', sep = '.'))
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
#' @param mir A required \code{character} vector of the microRNAs of interest.
#' These are the miRBase ID which are the official identifiers of the
#' widely used miRBase database, \url{http://www.mirbase.org/}.
#' @inheritParams stat_make
#' @inheritParams stat_collect
#'
#' @return A tidy \code{data.frame} of four columns. \code{mirna_base} is the
#' microRNA miRBase IDs, \code{feature} is the features/genes, \code{cor}
#' is the corresponding expression correlations and \code{study} is TCGA
#' study ID.
#' 
#' @examples 
#' # locate the testset file and connect
#' fl <- system.file('extdata', 'cRegulome.db', package = 'cRegulome')
#' conn <- RSQLite::dbConnect(RSQLite::SQLite(), fl)
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
#' @importFrom DBI dbListFields
#' 
#' @export
get_mir <- function(conn, mir, study, min_abs_cor, max_num,
                    targets_only = FALSE) {
    if(missing(mir)) {
        stop("User should provide at least one microRNA ID")
    } else if(!is.character(mir)) {
        stop("mir should be a character vector")
    } else {
        tf_id <- as.character(mir)
    }
    
    if(!missing(study)) {
        if(!is.character(study)) {
            stop("Study should be a character vector.")
        }
        study <- as.character(study)
    } else {
        study <- dbListFields(conn, 'cor_mir')[-1:-2]
    }
    
    if(!missing(min_abs_cor)) {
        if(!is.numeric(min_abs_cor)) {
            stop("min_abs_cor should be a numeric between 0 and 1.")
            }
        if(min_abs_cor > 1 || min_abs_cor < 0) {
            stop('min_abs_cor should be a numeric between 0 and 1.')
        }
    }
    if(!missing(max_num)) {
        if(!is.numeric(max_num)) {
            stop('max_num should be a positive integer.')
        }
        if(max_num < 1) {
            stop('max_num should be a positive integer.')
        }
    }
    
    # construct and excute query
    ll1 <- list()
    ll2 <- list()
    
    for(m in 1:length(mir)) {
        for(s in 1:length(study)) {
             stat <- stat_make(mir[m],
                              study = study[s],
                              min_abs_cor = min_abs_cor,
                              max_num = max_num)
             df <- stat_collect(conn,
                                study = study[s],
                                stat)
             ll2[[s]] <- df
        }
        ll1[[m]] <- ll2
    }
    
    # make data.frame
    ll <- unlist(ll1, recursive = FALSE)
    dat <- do.call('rbind', ll)
    
    # return cor to the -1:1 range
    dat$cor <- dat$cor/100
    
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
#' @inheritParams stat_make
#' @inheritParams stat_collect
#'
#' @return A tidy \code{data.frame} of four columns. \code{tf} is the official
#' gene symbols of the genes contains the transcription factor, \code{feature}
#' is the features/genes, cor is the corresponding expression correlations
#' and \code{study} is TCGA study ID.
#' 
#' @examples
#' # locate the testset file and connect
#' fl <- system.file('extdata', 'cRegulome.db', package = 'cRegulome')
#' conn <- RSQLite::dbConnect(RSQLite::SQLite(), fl)
#' 
#' \dontrun{
#' # get transcription factors correlations in all studies
#' get_tf(conn,
#'         tf = 'LEF1')
#' }
#' 
#' # get correlations in a particular study
#' get_tf(conn,
#'        tf = 'LEF1',
#'        study = '"STES*"')
#' 
#' # enter a custom query with different arguments
#' get_tf(conn,
#'        tf = 'LEF1',
#'        study = '"STES*"',
#'        min_abs_cor = .3,
#'        max_num = 5)
#' 
#' @importFrom DBI dbListFields
#' 
#' @export
get_tf <- function(conn, tf, study, min_abs_cor, max_num,
                    targets_only = FALSE) {
    if(missing(tf)) {
        stop("User should provide at least one TF ID")
    } else if(!is.character(tf)) {
        stop("tf should be a character vector")
    } else {
        tf_id <- as.character(tf)
    }
    
    if(!missing(study)) {
        if(!is.character(study)) {
            stop("Study should be a character vector.")
        }
        study <- as.character(study)
    } else {
        study <- dbListFields(conn, 'cor_tf')[-1:-2]
    }
    
    if(!missing(min_abs_cor)) {
        if(!is.numeric(min_abs_cor)) {
            stop("min_abs_cor should be a numeric between 0 and 1.")
        }
        if(min_abs_cor > 1 || min_abs_cor < 0) {
            stop('min_abs_cor should be a numeric between 0 and 1.')
        }
    }
    
    if(!missing(max_num)) {
        if(!is.numeric(max_num)) {
            stop('max_num should be a positive integer.')
        }
        if(max_num < 1) {
            stop('max_num should be a positive integer.')
        }
    }
    
    # construct and excute query
    ll1 <- list()
    ll2 <- list()
    
    for(m in 1:length(tf)) {
        for(s in 1:length(study)) {
            stat <- stat_make(tf[m],
                              study = study[s],
                              min_abs_cor = min_abs_cor,
                              max_num = max_num,
                              type = 'tf')
            df <- stat_collect(conn,
                               study = study[s],
                               stat,
                               type = 'tf')
                ll2[[s]] <- df
        }
        ll1[[m]] <- ll2
    }
    
    # make data.frame
    ll <- unlist(ll1, recursive = FALSE)
    dat <- do.call('rbind', ll)
    
    # return cor to the -1:1 range
    dat$cor <- dat$cor/100
    
    # return dat
    return(dat)
}

#' Make A SQL statment
#' 
#' Not meant to be called direclty by the user.
#'
#' @param reg A \code{character} vector of one or more regulator ID.
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
#' @param type A \code{character} string. Either 'mir' of 'tf'. Used to define
#' columns and tables names.
#' 
#' @examples 
#' stat_make(reg = 'hsa-let-7g',
#'           study = 'STES')
#'           
#' stat_make(reg = 'hsa-let-7g',
#'           study = 'STES',
#'           min_abs_cor = .3)
#'           
#' stat_make(reg = 'hsa-let-7g',
#'           study = 'STES',
#'           min_abs_cor = .3,
#'           max_num = 5)
#'           
#' @return A character string
#' 
#' @export
stat_make <- function(reg, study, min_abs_cor, max_num, targets_only = FALSE,
                      type = 'mir') {
    # define columns and tables based on type
    # column name
    id <- switch(type,
                 'mir' = 'mirna_base',
                 'tf' = 'tf'
    )
    
    # correlation table
    cor_tab <- switch(type,
                      'mir' = 'cor_mir',
                      'tf' = 'cor_tf' 
    )
    
    # targets table
    targets_tab <- switch(type,
                          'mir' = 'targets_mir',
                          'tf' = 'targets_tf')
    
    # make statement
    ## main select statement
    main <- paste0(
        'SELECT ',
        cor_tab, '.', id, ', ',
        cor_tab, '.feature, ',
        cor_tab, '.', study,
        ' FROM ', cor_tab
    )
    
    ## select targets only
    if(targets_only) {
        join <- paste0(
            'INNER JOIN ',
            targets_tab, ' ON ',
            cor_tab, '.', id, ' = ', targets_tab, '.', id,
            ' AND ',
            cor_tab, '.feature', ' = ', targets_tab, '.feature'  
        )
        main <- paste(main, join)
    }
    
    ## filter one regulator
    whr <- paste0(
        'WHERE ', cor_tab, '.', id, '=', '"', reg, '"'
    )
    main <- paste(main, whr)
    
    ## minimum value
    if(!missing(min_abs_cor)) {
        fltr1 <- paste0(
            'AND ', 'ABS(', cor_tab, '.', study, ')', ' > ', abs(min_abs_cor) * 100
        )
        main <- paste(main, fltr1)
    }
    
    ## limit returned entries 
    if(!missing(max_num)) {
        fltr2 <- paste0(
            'ORDER BY ABS(', cor_tab, '.', study,') DESC ', 
            'LIMIT ', max_num
        )
        main <- paste(main, fltr2)
    }
    
    # return statement
    return(main)
}

#' Collect data from SQLite database
#' 
#' Not meant to be called direclty by the user.
#'
#' @param conn A connection such as this returned by 
#' \code{\link[DBI]{dbConnect}}
#' @inheritParams stat_make
#' @param stat A string such as this returned by \code{\link{stat_make}}
#'
#' @return A data.frame
#' 
#' @importFrom RSQLite dbGetQuery
#' 
#' @export
stat_collect <- function(conn, study, stat, type = 'mir') {
    # define colum name based on type
    id <- switch (type,
        'mir' = 'mirna_base',
        'tf' = 'tf'
    )
    
    # get query
    df <- dbGetQuery(conn, stat)
    
    if(nrow(df) > 0) {
        # add study column
        df$study <- study
        
        # rename columns
        names(df) <- c(id, 'feature', 'cor', 'study')
    } 
    
    # return data.frame
    return(df)
}
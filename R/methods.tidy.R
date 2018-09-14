#' Tidy \link{cmicroRNA} and \link{cTF} objects
#'
#' @inheritParams cor_plot
#'
#' @return A tidy \code{data.frame} of four columns. \code{mirna_base} or
#' \code{tf}is the microRNA miRBase IDs, \code{feature} is the features/genes,
#' \code{cor} is the corresponding expression correlations and \code{study}
#' is TCGA study ID.
#' 
#' @examples 
#' # locate the testset file and connect
#' fl <- system.file('extdata', 'cRegulome.db', package = 'cRegulome')
#' conn <- RSQLite::dbConnect(RSQLite::SQLite(), fl)
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
#' # convert cmicroRNA object to a tidy data.frame
#' tidy_cmir <- cor_tidy(cmir)
#'
#' @export
cor_tidy <- function(ob) {
    UseMethod('cor_tidy')
}

#' @export
cor_tidy.cmicroRNA <- function(ob) {
    # loop over studies and microRNAs
    # convert matrix to data.frame
    ll <- list()
    for(i in length(ob$corr)) {
        mat <- ob$corr[[i]]
        mir <- colnames(mat)
        dfs <- list()
        for(t in 1:length(mir)) {
            df <- data.frame(
                cor = mat[, mir[t]],
                mirna_base = mir[t],
                feature = rownames(mat),
                stringsAsFactors = FALSE
            )
            
            dfs[[t]] <- df
        }
        dfs <- do.call('rbind', dfs)
        dfs$study <- ob$studies[i]
        ll[[i]] <- dfs
    }
    
    # bind all data.frames into a single one
    ll <- do.call('rbind', ll)
    
    # order the columns and remove row.names
    dat <- ll[, c('mirna_base', 'feature', 'cor', 'study')]
    row.names(dat) <- NULL
    dat <- dat[!is.na(dat$cor),]
    
    # return data.frame
    return(dat)
}

#' @export
cor_tidy.cTF <- function(ob) {
    # loop over studies and tfs 
    # convert matrix to data.frame
    ll <- list()
    for(i in length(ob$corr)) {
        mat <- ob$corr[[i]]
        tf <- colnames(mat)
        dfs <- list()
        for(t in 1:length(tf)) {
            df <- data.frame(
                cor = mat[, tf[t]],
                tf = tf[t],
                feature = rownames(mat),
                stringsAsFactors = FALSE
            )
            
            dfs[[t]] <- df
        }
        dfs <- do.call('rbind', dfs)
        dfs$study <- ob$studies[i]
        ll[[i]] <- dfs
    }
    
    # bind all data.frames into a single one
    ll <- do.call('rbind', ll)
    
    # order the columns and remove row.names
    dat <- ll[, c('tf', 'feature', 'cor', 'study')]
    row.names(dat) <- NULL
    dat <- dat[!is.na(dat$cor),]
    
    # return data.frame
    return(dat)
}

#' Make an igraph object
#' 
#' An \code{igraph} object of from \link{cmicroRNA} or \link{cTF} 
#' objects. 
#' 
#' @inheritParams cor_plot
#'
#' @return An \code{igraph} object 
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
#'                mir = c('hsa-let-7g', 'hsa-let-7i'),
#'                study = 'STES')
#' 
#' # make a cmicroRNA object   
#' cmir <- cmicroRNA(dat)
#' 
#' # print object
#' cor_igraph(cmir)
#' 
#' @importFrom igraph graph_from_data_frame
#' 
#' @export
cor_igraph <- function(ob) {
    UseMethod('cor_igraph')
}

#' @export
cor_igraph.cmicroRNA <- function(ob) {
    # get a tidy data.frame of the object
    dat <- cor_tidy(ob)
    
    # make edges
    edgs <- data.frame(
        from = dat$mirna_base,
        to = dat$feature,
        weight = abs(dat$cor),
        stringsAsFactors = FALSE
    )
    
    # make vertices
    vrtcs <- list(data.frame(
        id = unique(edgs$from),
        type = 'microRNA'
    ),
    data.frame(
        id = unique(edgs$to),
        type = 'gene'
    ))
    
    vrtcs <- do.call('rbind', vrtcs)
    
    # make graph
    g <- graph_from_data_frame(d = edgs,
                               directed = FALSE,
                               vrtcs)
    # return graph
    return(g)
}

#' @export
cor_igraph.cTF <- function(ob) {
    # get a tidy data.frame of the object
    dat <- cor_tidy(ob)
    
    # make edges
    edgs <- data.frame(
        from = dat$tf,
        to = dat$feature,
        weight = abs(dat$cor),
        stringsAsFactors = FALSE
    )
    
    # make vertices
    vrtcs <- list(data.frame(
        id = unique(edgs$from),
        type = 'TF'
    ),
    data.frame(
        id = unique(edgs$to),
        type = 'gene'
    ))
    
    vrtcs <- do.call('rbind', vrtcs)
    
    # make graph
    g <- graph_from_data_frame(d = edgs,
                               directed = FALSE,
                               vrtcs)
    
    # return graph
    return(g)
}

#' Prepare correlation data for plotting
#' 
#' Not meant to be called direclty by the user.
#'
#' @param ob A \link{cmicroRNA} or \link{cTF} object such as this returned by
#' calling \link{cmicroRNA} or \link{cTF}.
#' @param study A \code{character} vector of The Cancer Genome Atlas (TCGA)
#' study identifiers. To view the available studies in TCGA project,
#' \url{https://tcga-data.nci.nih.gov/docs/publications/tcga}. When left to
#' default \code{NULL} all available studies will be included.
#' @param add_dir A \code{logical} default TRUE for whether to add a column
#' called Direction that has the direction of the correlation; positive or 
#' negative.
#' @param add_corr A \code{logical} default TRUE for whether to add a colum
#' called Correlation that has the absolute value of the correlation
#'
#' @return A \code{data.frame}
#' 
#' @export
cor_prep <- function(ob, study, add_dir = TRUE, add_corr = TRUE) {
    # check the validity of the input study
    studies <- ob$studies
    
    # subset the data to the input study
    if(length(studies) > 1 & missing(study)) {
        warning('ob has multiple studies. First study will be used.')
        study <- studies[[1]]
    } else if (missing(study)){
        study <- studies[[1]]
    }
    
    # prepare data
    dat <- cor_tidy(ob)
    dat <- dat[dat$study == study,]
    dat <- dat[!is.na(dat$cor),]
    
    # create correlation and direction variables
    dat['Correlation'] <- abs(dat$cor)
    dat['Direction'] <- ifelse(dat$cor > 0, 'Positive', 'Negative')
    
    return(dat)
}

#' @export
print.cmicroRNA <- function(x, ...) {
    p <- paste(
        'A cmicroRNA object: microRNA-gene correlations in Cancer\n',
        'Contains:\n',
        length(x$studies), 'Cancer study/ies:', paste(x$studies,
                                                      collapse = ' '),
        '\n',
        length(x$microRNA), 'microRNA/s:', paste(x$microRNA,
                                                 collapse = ' '),
        '\n',
        length(x$features), 'features:', paste(x$features[1:5],
                                               collapse = ' '),
        '\n')
    cat(p)
}

#' @export
print.cTF <- function(x, ...) {
    p <- paste(
        'A cTF object: transcription factor-gene correlations in Cancer\n',
        'Contains:\n',
        length(x$studies), 'Cancer study/ies:', paste(x$studies,
                                                      collapse = ' '),
        '\n',
        length(x$TF), 'Transcription factor/s:', paste(x$TF,
                                                       collapse = ' '),
        '\n',
        length(x$features), 'features:', paste(x$features[1:5],
                                               collapse = ' '),
        '\n')
    cat(p)
}

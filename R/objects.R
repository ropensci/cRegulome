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
#' @export
cmicroRNA <- function(dat_mir){
    # extract items of the list from the data.frame
    microRNA <- unique(dat_mir$mirna_base)
    features <- unique(dat_mir$feature)
    studies <- unique(dat_mir$study)

    # reshape the data.frame/s 
    # microRNA in columns and feature in rows
    dfs <- split(dat_mir, dat_mir$study)
    
    corr <- list()
    for(l in 1:length(studies)) {
        df <- dfs[[l]]
        
        ll <- list()
        for(i in 1:length(microRNA)) {
            d <- df[df$mirna_base == microRNA[i],]
            m <- matrix(d$cor, ncol = 1)
            rownames(m) <- d$feature
            colnames(m) <- microRNA[i]
            ll[[i]] <- m
        }
        
        mat <- Reduce(function(x, y) {
            m <- merge(x, y, by = 'row.names', all = TRUE)
            if(colnames(m)[1] == 'Row.names') {
                rownames(m) <- m[, 1]
                m <- m[, -1]
                m
            }
        }, ll)
        
        corr[[l]] <- mat
    }
    
    names(corr) <- studies
    
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
#' # locate the testset file and connect
#' fl <- system.file('extdata', 'cRegulome.db', package = 'cRegulome')
#' conn <- RSQLite::dbConnect(RSQLite::SQLite(), fl)
#' 
#' # enter a custom query with different arguments
#' dat <- get_tf(conn,
#'               tf = 'LEF1',
#'               study = 'STES',
#'               min_abs_cor = .3,
#'               max_num = 5)
#' 
#' # make a cTF object   
#' ctf <- cTF(dat)
#' 
#' @export
cTF <- function(dat_tf){
    # extract items of the list from the data.frame
    TF <- unique(dat_tf$tf)
    features <- unique(dat_tf$feature)
    studies <- unique(dat_tf$study)


        
    # reshape the data.frame/s 
    # tf in columns and feature in rows
    dfs <- split(dat_tf, dat_tf$study)
    
    corr <- list()
    for(l in 1:length(studies)) {
        df <- dfs[[l]]
        
        ll <- list()
        for(i in 1:length(TF)) {
            d <- df[df$tf == TF[i],]
            m <- matrix(d$cor, ncol = 1)
            rownames(m) <- d$feature
            colnames(m) <- TF[i]
            ll[[i]] <- m
        }
        
        mat <- Reduce(function(x, y) {
            m <- merge(x, y, by = 'row.names', all = TRUE)
            if(colnames(m)[1] == 'Row.names') {
                rownames(m) <- m[, 1]
                m <- m[, -1]
            m
            }
            }, ll)
        
        corr[[l]] <- mat
    }
    
    names(corr) <- studies
    
    # object structure and class name
    structure(list(
        TF = TF,
        features = features,
        studies = studies,
        corr = corr
    ),
    class = 'cTF')
}

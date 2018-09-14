#' Plot method for \link{cmicroRNA} and \link{cTF} objects
#'
#' A dot plot of microRNA/TF correlation in a single study of TCGA. When the
#' object \link{cmicroRNA}/\link{cTF} contains more than one TCGA studies, the
#' argument \code{study} is a requirement.
#' 
#' @inheritParams cor_prep
#' @param ... Other options
#'
#' @return A \code{ggplot} object of a dot plot of the correlation values 
#' between genes and microRNAs or transcription factors in a TCGA study.
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
#' # print object
#' cor_plot(cmir)
#' 
#' @import ggplot2
#' 
#' @export
cor_plot <- function(ob, study, ...) {
    UseMethod('cor_plot')
}

#' @export
cor_plot.cmicroRNA <- function(ob, study, ...) {
    # prepare data for plotting
    dat <- cor_prep(ob, study = study)
    
    # plot 
    gg <- ggplot(dat, aes_string(x = 'mirna_base',
                                 y = 'feature',
                                 size = 'Correlation',
                                 color = 'Direction')) +
        geom_point() +
        theme_light()
    
    # return plot
    return(gg)
}

#' @export
cor_plot.cTF <- function(ob, study, ...) {
    # prepare data for plotting
    dat <- cor_prep(ob, study = study)
    
    # plot 
    gg <- ggplot(dat, aes_string(x = 'tf',
                                 y = 'feature',
                                 size = 'Correlation',
                                 color = 'Direction')) +
        geom_point() +
        theme_light()
    
    # return plot
    return(gg)
}

#' Venn Diagram of microRNA or transcription factor correlated features
#'
#' Count and plot the numbers of microRNA correlated features in
#' \code{cmicroRNA} object.
#'
#' @inheritParams cor_plot
#' @return A venn diagram with a circle or an ellipses for each microRNA and
#' the number of correlated features.
#' 
#' @examples 
#' # locate the testset file and connect
#' fl <- system.file('extdata', 'cRegulome.db', package = 'cRegulome')
#' conn <- RSQLite::dbConnect(RSQLite::SQLite(), fl)
#' 
#' # enter a custom query with different arguments
#' dat <- get_mir(conn,
#'                mir = c('hsa-let-7g', 'hsa-let-7i'),
#'                study = 'STES')
#' 
#' # make a cmicroRNA object   
#' cmir <- cmicroRNA(dat)
#' 
#' # make graph
#' cor_venn_diagram(cmir)
#' 
#' @importFrom VennDiagram venn.diagram
#' @importFrom grid grid.draw
#' 
#' @export
cor_venn_diagram <- function(ob, study, ...) {
    UseMethod('cor_venn_diagram')
}

#' @export
cor_venn_diagram.cmicroRNA <- function(ob, study, ...) {
    # stop if one set of microRNAs
    mir <- ob$microRNA
    if(length(mir) == 1) {
        stop(
            paste('ob contains data for one microRNA.',
                'Venn diagram is used to compare at least two sets.'
            )
        )
    }
    
    # prepare data for plotting
    dat <- cor_prep(ob, study = study, add_dir = FALSE, add_corr = FALSE)
    
    # make a named list of features
    dat <- with(dat, split(feature, mirna_base))
    
    # generate plot
    pp <- venn.diagram(dat,
                       imagetype = 'png', 
                       filename = NULL, ...)
    
    grid.draw(pp)
}

#' @export
cor_venn_diagram.cTF <- function(ob, study, ...) {
    # stop if one set of TF
    tf <- ob$TF
    if(length(tf) == 1) {
        stop(
            paste('ob contains data for one TF.',
                  'Venn diagram is used to compare at least two sets.'
            )
        )
    }
    
    # prepare data for plotting
    dat <- cor_prep(ob, study = study, add_dir = FALSE, add_corr = FALSE)
    
    
    # make a named list of features
    dat <- with(dat, split(feature, tf))
    
    # generate plot
    pp <- venn.diagram(dat,
                       imagetype = 'png', 
                       filename = NULL, ...)
    
    grid.draw(pp)
}

#' \code{\link[UpSetR]{upset}} plot of microRNA or tf sets
#'
#' \code{\link[UpSetR]{upset}} of sets of microRNAs or transcription
#' factors and their correlated features in a TCGA study.
#'
#' @inheritParams cor_plot
#'
#' @return An \code{\link[UpSetR]{upset}} plot
#' 
#' @examples  
#' # locate the testset file and connect
#' fl <- system.file('extdata', 'cRegulome.db', package = 'cRegulome')
#' conn <- RSQLite::dbConnect(RSQLite::SQLite(), fl)
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
#' cor_upset(cmir)
#' 
#' @importFrom UpSetR upset
#' 
#' @export
cor_upset <- function(ob, study, ...) {
    UseMethod('cor_upset')
}

#' @export
cor_upset.cmicroRNA <- function(ob, study, ...) {
    # stop if one set of microRNAs
    mir <- ob$microRNA
    if(length(mir) == 1) {
        stop(
            paste('ob contains data for one microRNA.',
                  'upset plot is used to compare at least two sets.'
            )
        )
    }
    # check the validity of the input study
    studies <- ob$studies
    
    # subset the data to the input study
    if(length(studies) > 1 & missing(study)) {
        warning('ob has multiple studies. First study will be used.')
        study <- studies[[1]]
    } else if (missing(study)){
        study <- studies[[1]]
    }
    
    # reshape data
    dat <- ob$corr[[study]]
    dat <- apply(dat, 2, function(x) ifelse(is.na(x), 0, 1))
    dat <- cbind(feature = rownames(dat), as.data.frame(dat))
    rownames(dat) <- NULL
    
    # generate plot
    upset(dat, ...)
}

#' @export
cor_upset.cTF <- function(ob, study, ...) {
    # stop if one set of TF
    tf <- ob$TF
    if(length(tf) == 1) {
        stop(
            paste('ob contains data for one TF.',
                  'upset plot is used to compare at least two sets.'
            )
        )
    }
    
    # check the validity of the input study
    studies <- ob$studies
    
    # subset the data to the input study
    if(length(studies) > 1 & missing(study)) {
        warning('ob has multiple studies. First study will be used.')
        study <- studies[[1]]
    } else if (missing(study)){
        study <- studies[[1]]
    }
    
    # reshape data
    dat <- ob$corr[[study]]
    dat <- apply(dat, 2, function(x) ifelse(is.na(x), 0, 1))
    dat <- cbind(feature = rownames(dat), as.data.frame(dat))
    rownames(dat) <- NULL
    
    # generate plot
    upset(dat, ...)
}

#' A histogram of the correlations of microRNA or tf sets
#'
#' Plot a \code{\link[graphics]{hist}} of sets of microRNAs or transcription
#' factors-gene correlations in a TCGA study.
#'
#' @inheritParams cor_plot
#'
#' @return An \code{\link[graphics]{hist}} plot of the correlations values 
#' between genes a microRNA or a transcription factor in a TCGA study
#' 
#' @examples 
#' # locate the testset file and connect
#' fl <- system.file('extdata', 'cRegulome.db', package = 'cRegulome')
#' conn <- RSQLite::dbConnect(RSQLite::SQLite(), fl)
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
#' cor_hist(cmir)
#' 
#' @importFrom graphics hist
#' 
#' @export
cor_hist <- function(ob, study, ...) {
    UseMethod('cor_hist')
}

#' @export
cor_hist.cmicroRNA <- function(ob, study, ...) {
    # prepare data for plotting
    dat <- cor_prep(ob, study = study, add_dir = FALSE, add_corr = FALSE)
    
    # generate plot
    dat <- dat$cor
    hist(dat, ...)
}

#' @export
cor_hist.cTF <- function(ob, study, ...) {
    # prepare data for plotting
    dat <- cor_prep(ob, study = study, add_dir = FALSE, add_corr = FALSE)
    
    # generate plot
    dat <- dat$cor
    hist(dat, ...)
}

#' A joy plot of correlation of microRNA or tf sets
#'
#' A \code{\link{ggridges}} joy plot of sets of microRNAs or transcription
#' factors-gene correlations in a TCGA study.
#'
#' @inheritParams cor_plot
#'
#' @return An \code{\link{ggridges}} plot object
#' 
#' @examples 
#' # locate the testset file and connect
#' fl <- system.file('extdata', 'cRegulome.db', package = 'cRegulome')
#' conn <- RSQLite::dbConnect(RSQLite::SQLite(), fl)
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
#' cor_joy(cmir)
#' 
#' @importFrom ggridges geom_density_ridges
#' @import ggplot2 
#' 
#' @export
cor_joy <- function(ob, study, ...) {
    UseMethod('cor_joy')
}

#' @export
cor_joy.cmicroRNA <- function(ob, study, ...) {
    # prepare data for plotting
    dat <- cor_prep(ob, study = study, add_dir = FALSE, add_corr = FALSE)
    
    # generate plot
    gg <- ggplot(dat, aes_string(x = 'cor',
                          y = 'mirna_base')) +
        geom_density_ridges() +
        theme_light()
    
    return(gg)
}

#' @export
cor_joy.cTF <- function(ob, study, ...) {
    # prepare data for plotting
    dat <- cor_prep(ob, study = study, add_dir = FALSE, add_corr = FALSE)
    
    # generate plot
    gg <- ggplot(dat, aes_string(x = 'cor',
                          y = 'tf')) +
        geom_density_ridges() +
        theme_light()
    
    return(gg)
}

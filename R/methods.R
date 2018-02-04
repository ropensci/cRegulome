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
        length(x$tf), 'Transcription factor/s:', paste(x$TF,
                                                       collapse = ' '),
        '\n',
        length(x$features), 'features:', paste(x$features[1:5],
                                               collapse = ' '),
        '\n')
    cat(p)
}

#' Plot method for \link{cmicroRNA} and \link{cTF} objects
#'
#' A dot plot of microRNA/TF correlation in a single study of TCGA. When the
#' object \link{cmicroRNA}/\link{cTF} contains more than one TCGA studies, the
#' argument \code{study} is a requirement.
#'
#' @param ob A \link{cmicroRNA} or \link{cTF} object such as this returned by
#' calling \link{cmicroRNA} or \link{cTF}.
#' @param study A \code{character} vector of The Cancer Genome Atlas (TCGA)
#' study identifiers. To view the available studies in TCGA project,
#' \url{https://tcga-data.nci.nih.gov/docs/publications/tcga}. When left to
#' default \code{NULL} all available studies will be included.
#' @param ... Other options
#'
#' @return A \code{ggplot} object of a dot plot of the correlation values 
#' between genes and microRNAs or transcription factors in a TCGA study.
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
#' # print object
#' cor_plot(cmir)
#' 
#' @importFrom magrittr %>%
#' 
#' @export
cor_plot <- function(ob, study = NULL, ...) {
    UseMethod('cor_plot')
}

#' @export
cor_plot.cmicroRNA <- function(ob, study = NULL, ...) {
    # check the validity of the input study
    if(length(ob$studies) > 1) {
        if(is.null(study) || length(study) > 1) {
            stop('User should provide a single study to plot.')
        }
    }
    
    # subset the data to the input study
    if(is.data.frame(ob$corr)) {
        # object contains a single study
        dat <- ob$corr
    } else {
        # object contains multiple studies
        dat <- ob$corr[[study]]
    }
    
    
    
    # convet back to tidy format
    dat <- dat %>%
        tidyr::gather(mirna_base, cor, -feature) %>%
        stats::na.omit()
    
    # create correlation and direction variables
    # plot 
    gg <- dat %>%
        dplyr::mutate(Correlation = abs(cor),
                      Direction = ifelse(cor > 0, 'Positive', 'Negative')) %>%
        ggplot2::ggplot(ggplot2::aes_string(x = 'mirna_base',
                                            y = 'feature',
                                            size = 'Correlation',
                                            color = 'Direction')) +
        ggplot2::geom_point() +
        ggplot2::theme_light()
    
    return(gg)
}

#' @export
cor_plot.cTF <- function(ob, study = NULL, ...) {
    # check the validity of the input study
    if(length(ob$studies) > 1) {
        if(is.null(study) || length(study) > 1) {
            stop('User should provide a single study to plot.')
        }
    }
    
    # subset the data to the input study
    if(is.data.frame(ob$corr)) {
        # object contains a single study
        dat <- ob$corr
    } else {
        # object contains multiple studies
        dat <- ob$corr[[study]]
    }
    
    
    
    # convet back to tidy format
    dat <- dat %>%
        tidyr::gather(tf, cor, -feature) %>%
        stats::na.omit()
    
    # create correlation and direction variables
    # plot 
    gg <- dat %>%
        dplyr::mutate(Correlation = abs(cor),
                      Direction = ifelse(cor > 0, 'Positive', 'Negative')) %>%
        ggplot2::ggplot(ggplot2::aes_string(x = 'tf',
                                            y = 'feature',
                                            size = 'Correlation',
                                            color = 'Direction')) +
        ggplot2::geom_point() +
        ggplot2::theme_light()
    
    
    return(gg)
}

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
#' # convert cmicroRNA object to a tidy data.frame
#' tidy_cmir <- cor_tidy(cmir)
#' 
#' @importFrom magrittr %>%
#'
#' @export
cor_tidy <- function(ob) {
    UseMethod('cor_tidy')
}

#' @export
cor_tidy.cmicroRNA <- function(ob) {
    
    
    # convet back to tidy format
    if(length(ob$studies) == 1) {
        # object contains a single study
        dat <- ob$corr %>%
            tidyr::gather(mirna_base, cor, -feature) %>%
            stats::na.omit() %>%
            dplyr::mutate(study = ob$studies)
    } else {
        # object with multiple studies
        dat <- purrr::map(ob$corr, function(x) {
            x %>%
                tidyr::gather(mirna_base, cor, -feature) %>%
                stats::na.omit()
        }) %>%
            dplyr::bind_rows(.id = 'study') %>%
            dplyr::select(2:4, study)
    }
    
    return(dat)
}

#' @export
cor_tidy.cTF <- function(ob) {
    
    
    # convet back to tidy format
    if(length(ob$studies) == 1) {
        # object contains a single study
        dat <- ob$corr %>%
            tidyr::gather(tf, cor, -feature) %>%
            stats::na.omit() %>%
            dplyr::mutate(study = ob$studies)
    } else {
        # object with multiple studies
        dat <- purrr::map(ob$corr, function(x) {
            x %>%
                tidyr::gather(tf, cor, -feature) %>%
                stats::na.omit()
        }) %>%
            dplyr::bind_rows(.id = 'study') %>%
            dplyr::select(2:4, study)
    }
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
#' cor_venn_diagram(cmir)
#' 
#' @importFrom magrittr %>%
#' 
#' @export
cor_venn_diagram <- function(ob, study = NULL, ...) {
    UseMethod('cor_venn_diagram')
}

#' @export
cor_venn_diagram.cmicroRNA <- function(ob, study = NULL, ...) {
    # check the validity of the input study
    if(length(ob$studies) > 1) {
        if(is.null(study) || length(study) > 1) {
            stop('User should provide a single study to plot.')
        }
    }
    
    # subset the data to the input study
    if(is.data.frame(ob$corr)) {
        # object contains a single study
        dat <- ob$corr
    } else {
        # object contains multiple studies
        dat <- ob$corr[[study]]
    }
    
    
    
    # convert back to tidy format
    dat <- dat %>%
        tidyr::gather(mirna_base, cor, -feature) %>%
        stats::na.omit()
    
    # make a named list of features
    dat <- with(dat, split(feature, mirna_base))
    
    # generate plot
    pp <- VennDiagram::venn.diagram(dat,
                                    imagetype = 'png', 
                                    filename = NULL, ...)
    grid::grid.draw(pp)
}

#' @export
cor_venn_diagram.cTF <- function(ob, study = NULL, ...) {
    # check the validity of the input study
    if(length(ob$studies) > 1) {
        if(is.null(study) || length(study) > 1) {
            stop('User should provide a single study to plot.')
        }
    }
    
    # subset the data to the input study
    if(is.data.frame(ob$corr)) {
        # object contains a single study
        dat <- ob$corr
    } else {
        # object contains multiple studies
        dat <- ob$corr[[study]]
    }
    
    
    
    # convert back to tidy format
    dat <- dat %>%
        tidyr::gather(tf, cor, -feature) %>%
        stats::na.omit()
    
    # make a named list of features
    dat <- with(dat, split(feature, tf))
    
    # generate plot
    pp <- VennDiagram::venn.diagram(dat,
                                    imagetype = 'png', 
                                    filename = NULL, ...)
    grid::grid.draw(pp)
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
#' cor_upset(cmir)
#' 
#' @importFrom magrittr %>%
#' 
#' @export
cor_upset <- function(ob, study = NULL, ...) {
    UseMethod('cor_upset')
}

#' @export
cor_upset.cmicroRNA <- function(ob, study = NULL, ...) {
    # check the validity of the input study
    if(length(ob$studies) > 1) {
        if(is.null(study) || length(study) > 1) {
            stop('User should provide a single study to plot.')
        }
    }
    
    # subset the data to the input study
    if(is.data.frame(ob$corr)) {
        # object contains a single study
        dat <- ob$corr
    } else {
        # object contains multiple studies
        dat <- ob$corr[[study]]
    }
    
    
    
    # make a binary data.frame
    dat <- dat %>%
        dplyr::mutate_at(dplyr::vars(2:ncol(dat)),
                         function(x) x <- ifelse(is.na(x), 0, 1))
    
    # generate plot
    UpSetR::upset(dat, ...)
}

#' @export
cor_upset.cTF <- function(ob, study = NULL, ...) {
    # check the validity of the input study
    if(length(ob$studies) > 1) {
        if(is.null(study) || length(study) > 1) {
            stop('User should provide a single study to plot.')
        }
    }
    
    # subset the data to the input study
    if(is.data.frame(ob$corr)) {
        # object contains a single study
        dat <- ob$corr
    } else {
        # object contains multiple studies
        dat <- ob$corr[[study]]
    }
    
    
    
    # make a binary data.frame
    dat <- dat %>%
        dplyr::mutate_at(dplyr::vars(2:ncol(dat)),
                         function(x) x <- ifelse(is.na(x), 0, 1))
    
    # generate plot
    UpSetR::upset(dat, ...)
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
#' cor_hist(cmir)
#' 
#' @export
cor_hist <- function(ob, study = NULL, ...) {
    UseMethod('cor_hist')
}

#' @export
cor_hist.cmicroRNA <- function(ob, study = NULL, ...) {
    # check the validity of the input study
    if(length(ob$studies) > 1) {
        if(is.null(study) || length(study) > 1) {
            stop('User should provide a single study to plot.')
        }
    }
    
    # subset the data to the input study
    if(is.data.frame(ob$corr)) {
        # object contains a single study
        dat <- ob$corr
    } else {
        # object contains multiple studies
        dat <- ob$corr[[study]]
    }
    
    # generate plot
    dat <- unlist(dat[, -1])
    graphics::hist(dat, ...)
}

#' @export
cor_hist.cTF <- function(ob, study = NULL, ...) {
    # check the validity of the input study
    if(length(ob$studies) > 1) {
        if(is.null(study) || length(study) > 1) {
            stop('User should provide a single study to plot.')
        }
    }
    
    # subset the data to the input study
    if(is.data.frame(ob$corr)) {
        # object contains a single study
        dat <- ob$corr
    } else {
        # object contains multiple studies
        dat <- ob$corr[[study]]
    }
    
    # generate plot
    dat <- unlist(dat[, -1])
    graphics::hist(dat, ...)
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
#' cor_joy(cmir)
#' 
#' @importFrom magrittr %>%
#' 
#' @export
cor_joy <- function(ob, study = NULL, ...) {
    UseMethod('cor_joy')
}

#' @export
cor_joy.cmicroRNA <- function(ob, study = NULL, ...) {
    # check the validity of the input study
    if(length(ob$studies) > 1) {
        if(is.null(study) || length(study) > 1) {
            stop('User should provide a single study to plot.')
        }
    }
    
    # subset the data to the input study
    if(is.data.frame(ob$corr)) {
        # object contains a single study
        dat <- ob$corr
    } else {
        # object contains multiple studies
        dat <- ob$corr[[study]]
    }
    
    
    
    # convert back to tidy format
    dat <- dat %>%
        tidyr::gather(mirna_base, cor, -feature) %>%
        stats::na.omit()
    
    # generate plot
    gg <- dat %>%
        ggplot2::ggplot(ggplot2::aes_string(x = 'cor',
                                            y = 'mirna_base')) +
        ggridges::geom_density_ridges() +
        ggplot2::theme_light()
    
    return(gg)
}

#' @export
cor_joy.cTF <- function(ob, study = NULL, ...) {
    # check the validity of the input study
    if(length(ob$studies) > 1) {
        if(is.null(study) || length(study) > 1) {
            stop('User should provide a single study to plot.')
        }
    }
    
    # subset the data to the input study
    if(is.data.frame(ob$corr)) {
        # object contains a single study
        dat <- ob$corr
    } else {
        # object contains multiple studies
        dat <- ob$corr[[study]]
    }
    
    
    
    # convert back to tidy format
    dat <- dat %>%
        tidyr::gather(tf, cor, -feature) %>%
        stats::na.omit()
    
    # generate plot
    gg <- dat %>%
        ggplot2::ggplot(ggplot2::aes_string(x = 'cor',
                                            y = 'tf')) +
        ggridges::geom_density_ridges() +
        ggplot2::theme_light()
    
    return(gg)
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
        weight = abs(dat$cor)
    )
    
    # make vertices
    vrtcs <- list(microRNA = unique(edgs$from),
                  gene = unique(edgs$to))
    vrtcs <- reshape2::melt(vrtcs)
    names(vrtcs) <- c('id', 'type')
    
    # make graph
    g <- igraph::graph_from_data_frame(d = edgs,
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
        weight = abs(dat$cor)
    )
    
    # make vertices
    vrtcs <- list(TF = unique(edgs$from),
                  gene = unique(edgs$to))
    vrtcs <- unique(reshape2::melt(vrtcs))
    names(vrtcs) <- c('id', 'type')
    
    # make graph
    g <- igraph::graph_from_data_frame(d = edgs,
                                       directed = FALSE,
                                       vrtcs)
    
    # return graph
    return(g)
}

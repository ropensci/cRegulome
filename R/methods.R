#' Print method for \link{cmicroRNA} and \link{cTF} objects
#'
#' @param ob A \link{cmicroRNA} or \link{cTF} object such as this returned by
#' calling \link{cmicroRNA} or \link{cTF}.
#' @param ... Other argument to \code{\link[base]{print}}.
#'
#' @return Printed text
#'
#' @examples
#' # connect to test database file
#' db_file <- system.file("extdata", "cRegulome.db", package = "cRegulome")
#' conn <- DBI::dbConnect(RSQLite::SQLite(), db_file)
#'
#' # get data for 2 microRNAs in the ACC study
#' dat <- get_mir(conn,
#'     mir = c('hsa-let-7b', 'hsa-mir-134'),
#'     study = c('ACC', 'BLCA'),
#'     min_cor = .5)
#' DBI::dbDisconnect(conn)
#'
#' # convert to cmicroRNA object
#' ob <- cmicroRNA(dat)
#' print(ob)
#'
#' @export
print <- function(ob, ...) {
    UseMethod('print')
}

#' @export
print.cmicroRNA <- function(ob, ...) {
    p <- paste(
        'A cmicroRNA object: microRNA-gene correlations in Cancer\n',
        'Contains:\n',
        length(ob$studies), 'Cancer study/ies:', paste(ob$studies,
                                                    collapse = ' '),
        '\n',
        length(ob$microRNA), 'microRNA/s:', paste(ob$microRNA,
                                                collapse = ' '),
        '\n',
        length(ob$features), 'features:', paste(ob$features[1:5],
                                                collapse = ' '),
        '\n')
    cat(p)
}

#' @export
print.cTF <- function(ob, ...) {
    p <- paste(
        'A cTF object: transcription factor-gene correlations in Cancer\n',
        'Contains:\n',
        length(ob$studies), 'Cancer study/ies:', paste(ob$studies,
                                                        collapse = ' '),
        '\n',
        length(ob$tf), 'Transcription factor/s:', paste(ob$TF,
                                                        collapse = ' '),
        '\n',
        length(ob$features), 'features:', paste(ob$features[1:5],
                                                collapse = ' '),
        '\n')
    cat(p)
}

#' Plot method for \link{cmicroRNA} and \link{cTF} objects
#'
#' A dot plot of microRNA/TF correlation in a single study of TCGA. When the
#' object \link{cmicroRNA}/\link{cTF} contains more than one TCGA studies, the
#' argument \code{study} is a requirment.
#'
#' @param ob A \link{cmicroRNA} or \link{cTF} object such as this returned by
#' calling \link{cmicroRNA} or \link{cTF}.
#' @param study A \code{character} vector of The Cancer Genome Atlase (TCGA)
#' study identifiers. To view the available studies in TCGA project,
#' \url{https://tcga-data.nci.nih.gov/docs/publications/tcga}. When left to
#' defult \code{NULL} all available studies will be included.
#' @param ... Other options
#'
#' @return A \code{ggplot} object of dot plot.
#'
#' @examples
#' # connect to test database file
#' db_file <- system.file("extdata", "cRegulome.db", package = "cRegulome")
#' conn <- DBI::dbConnect(RSQLite::SQLite(), db_file)
#'
#' # get data for 2 microRNAs in the ACC study
#' dat <- get_mir(conn,
#'     mir = c('hsa-let-7b', 'hsa-mir-134'),
#'     study = c('ACC', 'BLCA'),
#'     min_cor = .5)
#' DBI::dbDisconnect(conn)
#'
#' # convert to cmicroRNA object
#' ob <- cmicroRNA(dat)
#' plot(ob, study = 'ACC')
#'
#' @export
plot <- function(ob, study = NULL, ...) {
    UseMethod('plot')
}

#' @export
plot.cmicroRNA <- function(ob, study = NULL, ...) {
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
  
    `%>%` <- dplyr::`%>%`
    
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
plot.cTF <- function(ob, study = NULL, ...) {
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
    
    `%>%` <- dplyr::`%>%`

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
#' @inheritParams plot
#'
#' @return A tidy \code{data.frame} of four columns. \code{mirna_base} or
#' \code{tf}is the microRNA miRBase IDs, \code{feature} is the features/genes,
#' \code{cor} is the corresponding expression correaltions and \code{study}
#' is TCGA study ID.
#'
#' @examples
#' # connect to test database file
#' db_file <- system.file("extdata", "cRegulome.db", package = "cRegulome")
#' conn <- DBI::dbConnect(RSQLite::SQLite(), db_file)
#'
#' # get data for 2 microRNAs in the ACC study
#' dat <- get_mir(conn,
#'     mir = c('hsa-let-7b', 'hsa-mir-134'),
#'     study = c('ACC', 'BLCA'),
#'     min_cor = .5)
#' DBI::dbDisconnect(conn)
#'
#' # convert to cmicroRNA object
#' ob <- cmicroRNA(dat)
#' dat <- tidy(ob)
#' dat[1:5,]
#'
#' @export
tidy <- function(ob) {
UseMethod('tidy')
}

#' @export
tidy.cmicroRNA <- function(ob) {
    `%>%` <- dplyr::`%>%`
    
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
tidy.cTF <- function(ob) {
    `%>%` <- dplyr::`%>%`
    
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
#' @inheritParams plot
#' @return A venn diagram with a ciccle or an ellipses for each microRNA and
#' the number of correlated features.
#'
#' @examples
#' # connect to test database file
#' db_file <- system.file("extdata", "cRegulome.db", package = "cRegulome")
#' conn <- DBI::dbConnect(RSQLite::SQLite(), db_file)
#'
#' # get data for 2 microRNAs in the ACC study
#' dat <- get_mir(conn,
#'     mir = c('hsa-let-7b', 'hsa-mir-134'),
#'     study = c('ACC', 'BLCA'),
#'     min_cor = .5)
#' DBI::dbDisconnect(conn)
#'
#' # convert to cmicroRNA object and plot
#' ob <- cmicroRNA(dat)
#' venn.diagram(ob, study = 'ACC')
#'
#' @export
venn.diagram <- function(ob, study = NULL, ...) {
    UseMethod('venn.diagram')
}

#' @export
venn.diagram.cmicroRNA <- function(ob, study = NULL, ...) {
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
    
    `%>%` <- dplyr::`%>%`

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
venn.diagram.cTF <- function(ob, study = NULL, ...) {
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
    
    `%>%` <- dplyr::`%>%`
    
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
#' @inheritParams plot
#'
#' @return An \code{\link[UpSetR]{upset}} plot
#'
#' @examples
#' # connect to test database file
#' db_file <- system.file("extdata", "cRegulome.db", package = "cRegulome")
#' conn <- DBI::dbConnect(RSQLite::SQLite(), db_file)
#'
#' # get data for 2 microRNAs in the ACC study
#' dat <- get_mir(conn,
#'     mir = c('hsa-let-7b', 'hsa-mir-134'),
#'     study = c('ACC', 'BLCA'))
#' DBI::dbDisconnect(conn)
#'
#' # convert to cmicroRNA object and plot
#' ob <- cmicroRNA(dat)
#' upset(ob, study = 'ACC')
#'
#' @export
upset <- function(ob, study = NULL, ...) {
    UseMethod('upset')
}

#' @export
upset.cmicroRNA <- function(ob, study = NULL, ...) {
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
    
    `%>%` <- dplyr::`%>%`
    
    # make a binary data.frame
    dat <- dat %>%
        dplyr::mutate_at(dplyr::vars(2:ncol(dat)),
                        function(x) x <- ifelse(is.na(x), 0, 1))
    
    # generate plot
    UpSetR::upset(dat, ...)
}

#' @export
upset.cTF <- function(ob, study = NULL, ...) {
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
    
    `%>%` <- dplyr::`%>%`
    
    # make a binary data.frame
    dat <- dat %>%
        dplyr::mutate_at(dplyr::vars(2:ncol(dat)),
                        function(x) x <- ifelse(is.na(x), 0, 1))
    
    # generate plot
    UpSetR::upset(dat, ...)
}

#' \code{\link{hist}} plot a histogram of microRNA or tf sets
#'
#' \code{\link{hist}} of sets of microRNAs or transcription
#' factors and their correlated features in a TCGA study.
#'
#' @inheritParams plot
#'
#' @return An \code{\link{hist}} plot
#'
#' @examples
#' # connect to test database file
#' db_file <- system.file("extdata", "cRegulome.db", package = "cRegulome")
#' conn <- DBI::dbConnect(RSQLite::SQLite(), db_file)
#'
#' # get data for 2 microRNAs in the ACC study
#' dat <- get_mir(conn,
#'     mir = c('hsa-let-7b', 'hsa-mir-134'),
#'     study = c('ACC', 'BLCA'))
#' DBI::dbDisconnect(conn)
#'
#' # convert to cmicroRNA object and plot
#' ob <- cmicroRNA(dat)
#' hist(ob, study = 'ACC')
#'
#' @export
hist <- function(ob, study = NULL, ...) {
    UseMethod('hist')
}

#' @export
hist.cmicroRNA <- function(ob, study = NULL, ...) {
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
    hist(dat, ...)
}

#' @export
hist.cTF <- function(ob, study = NULL, ...) {
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
    hist(dat, ...)
}

#' \code{\link{ggjoy}} joy plot of microRNA or tf sets
#'
#' \code{\link{ggjoy}} joy plot of sets of microRNAs or transcription
#' factors and their correlated features in a TCGA study.
#'
#' @inheritParams plot
#'
#' @return An \code{\link{ggjoy}} plot object
#'
#' @examples
#' # connect to test database file
#' db_file <- system.file("extdata", "cRegulome.db", package = "cRegulome")
#' conn <- DBI::dbConnect(RSQLite::SQLite(), db_file)
#'
#' # get data for 2 microRNAs in the ACC study
#' dat <- get_mir(conn,
#'     mir = c('hsa-let-7b', 'hsa-mir-134'),
#'     study = c('ACC', 'BLCA'))
#' DBI::dbDisconnect(conn)
#'
#' # convert to cmicroRNA object and plot
#' ob <- cmicroRNA(dat)
#' joy(ob, study = 'ACC')
#'
#' @export
joy <- function(ob, study = NULL, ...) {
    UseMethod('joy')
}

#' @export
joy.cmicroRNA <- function(ob, study = NULL, ...) {
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
    
    `%>%` <- dplyr::`%>%`
    
    # convert back to tidy format
    dat <- dat %>%
        tidyr::gather(mirna_base, cor, -feature) %>%
        stats::na.omit()
    
    # generate plot
    gg <- dat %>%
        ggplot2::ggplot(ggplot2::aes_string(x = 'cor',
                                            y = 'mirna_base')) +
        ggjoy::geom_joy() +
        ggplot2::theme_light()

    return(gg)
}

#' @export
joy.cTF <- function(ob, study = NULL, ...) {
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
    
    `%>%` <- dplyr::`%>%`
    
    # convert back to tidy format
    dat <- dat %>%
        tidyr::gather(tf, cor, -feature) %>%
        stats::na.omit()
    
    # generate plot
    gg <- dat %>%
        ggplot2::ggplot(ggplot2::aes_string(x = 'cor',
                                            y = 'tf')) +
        ggjoy::geom_joy() +
        ggplot2::theme_light()

    return(gg)
}
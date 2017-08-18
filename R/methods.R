#' Print method for \link{cmicroRNA} and \link{cTF} objects
#'
#' @param ob A \link{cmicroRNA} or \link{cTF} object such as this returned by
#'    calling \link{cmicroRNA} or \link{cTF}.
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
    p <- paste('A cmicroRNA object: microRNA expression correlations
             in Cancer\n',
               'Contains:\n',
               length(ob$studies), 'Cancer study/ies:', paste(ob$studies[1:5],
                                                              collapse = ' '),
               '\n',
               length(ob$microRNA), 'microRNA/s:', paste(ob$microRNA[1:5],
                                                         collapse = ' '),
               '\n',
               length(ob$features), 'features:', paste(ob$features[1:5],
                                                       collapse = ' '),
               '\n')
    cat(p)
}

#' @export
print.cTF <- function(ob, ...) {
    p <- paste('A cTF object: transcription factor expression correlations
             in Cancer\n',
               'Contains:\n',
               length(ob$studies), 'Cancer study/ies:', paste(ob$studies[1:5],
                                                              collapse = ' '),
               '\n',
               length(ob$tf), 'Transcription factor/s:', paste(ob$TF[1:5],
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
#'    calling \link{cmicroRNA} or \link{cTF}.
#' @param study A \code{character} vector of The Cancer Genome Atlase (TCGA)
#'    study identifiers. To view the available studies in TCGA project,
#'    \url{https://tcga-data.nci.nih.gov/docs/publications/tcga}. When left to
#'    defult \code{NULL} all available studies will be included.
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
    if(length(ob$studies) > 1 && length(study) != 1) {
        stop('User shouls provide a singl study to plot.')
    }

    if(is.data.frame(ob$corr)) {
        dat <- ob$corr
    } else {
        dat <- ob$corr[[study]]
    }
    `%>%` <- dplyr::`%>%`
    dat <- dat %>%
        tidyr::gather(mirna_base, cor, -feature) %>%
        stats::na.omit()
    gg <- dat %>%
        ggplot2::ggplot(ggplot2::aes(x = mirna_base, y = feature, size = abs(cor))) +
        ggplot2::geom_point()
    return(gg)
}

#' @export
plot.cTF <- function(ob, study = NULL, ...) {
    if(length(ob$studies) > 1 && length(study) != 1) {
        stop('User shouls provide a singl study to plot.')
    }

    if(is.data.frame(ob$corr)) {
        dat <- ob$corr
    } else {
        dat <- ob$corr[[study]]
    }
    `%>%` <- dplyr::`%>%`
    dat <- dat %>%
        tidyr::gather(tf, cor, -feature) %>%
        stats::na.omit()
    gg <- dat %>%
        ggplot2::ggplot(ggplot2::aes(x = tf, y = feature, size = abs(cor))) +
        ggplot2::geom_point()
    return(gg)
}

#' Tidy \link{cmicroRNA} and \link{cTF} objects
#'
#' @inheritParams plot
#'
#' @return A tidy \code{data.frame} of four columns. \code{mirna_base} or
#'    \code{tf}is the microRNA miRBase IDs, \code{feature} is the
#'    features/genes, \code{cor} is the corresponding expression
#'    correaltions and \code{study} is TCGA study ID.
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
    if(length(ob$studies) == 1) {
        dat <- ob$corr %>%
            tidyr::gather(mirna_base, cor, -feature) %>%
            stats::na.omit() %>%
            dplyr::mutate(study = ob$studies)
    } else {
        dat <- purrr::map(ob$corr, function(x) {
            x %>%
                tidyr::gather(mirna_base, cor, -feature) %>%
                stats::na.omit()
        }) %>%
            dplyr::bind_rows(.id = 'study') %>%
            dplyr::select(2:3, study)
    }
    return(dat)
}

#' @export
tidy.cTF <- function(ob) {
    `%>%` <- dplyr::`%>%`
    if(length(ob$studies) == 1) {
        dat <- ob$corr %>%
            tidyr::gather(tf, cor, -feature) %>%
            stats::na.omit() %>%
            dplyr::mutate(study = ob$studies)
    } else {
        dat <- purrr::map(ob$corr, function(x) {
            x %>%
                tidyr::gather(tf, cor, -feature) %>%
                stats::na.omit()
        }) %>%
            dplyr::bind_rows(.id = 'study') %>%
            dplyr::select(2:3, study)
    }
}

#' Venn Diagram of microRNA or transcription factor correlated features
#'
#' Count and plot the numbers of microRNA correlated features in
#' \code{cmicroRNA} object.
#'
#' @inheritParams plot
#' @return A venn diagram with a ciccle or an ellipses for each microRNA and
#'    the number of correlated features.
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
#' venn.diagram(ob, study = 'ACC', filename = 'mir.tiff')
#'
#' @export
venn.diagram <- function(ob, study = NULL, ...) {
    UseMethod('venn.diagram')
}

#' @rdname venn.daigram
#' @export
venn.diagram.cmicroRNA <- function(ob, study = NULL, ...) {
    if(length(ob$studies) > 1 && length(study) != 1) {
        stop('User shouls provide a singl study to plot.')
    }

    if(is.data.frame(ob$corr)) {
        dat <- ob$corr
    } else {
        dat <- ob$corr[[study]]
    }
    `%>%` <- dplyr::`%>%`
    dat <- dat %>%
        tidyr::gather(mirna_base, cor, -feature) %>%
        stats::na.omit()
    dat <- with(dat, split(feature, mirna_base))
    VennDiagram::venn.diagram(dat, ...)
}

#' @export
venn.diagram.cTF <- function(ob, study = NULL, ...) {
    if(length(ob$studies) > 1 && length(study) != 1) {
        stop('User should provide a singl study to plot.')
    }

    if(is.data.frame(ob$corr)) {
        dat <- ob$corr
    } else {
        dat <- ob$corr[[study]]
    }
    `%>%` <- dplyr::`%>%`
    dat <- dat %>%
        tidyr::gather(tf, cor, -feature) %>%
        stats::na.omit()
    dat <- with(dat, split(feature, tf))
    VennDiagram::venn.diagram(dat, ...)
}

#' \code{\link[UpSetR]{upset}} plot of microRNA or tf sets
#'
#' \code{\link[UpSetR]{upset}} of sets of microRNAs and their correlated
#' features in a TCGA study.
#'
#' @inheritParams plot
#' @return An \code{\link[UpSetR]{upset}} plot
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
    if(length(ob$studies) > 1 && length(study) != 1) {
        stop('User shouls provide a singl study to plot.')
    }

    if(is.data.frame(ob$corr)) {
        dat <- ob$corr
    } else {
        dat <- ob$corr[[study]]
    }
    `%>%` <- dplyr::`%>%`
    dat <- dat %>%
        dplyr::mutate_at(dplyr::vars(2:ncol(dat)),
                  function(x) x <- ifelse(is.na(x), 0, 1))
    UpSetR::upset(dat, ...)
}

#' @export
upset.cTF <- function(ob, study = NULL, ...) {
    if(length(ob$studies) > 1 && length(study) != 1) {
        stop('User shouls provide a singl study to plot.')
    }

    if(is.data.frame(ob$corr)) {
        dat <- ob$corr
    } else {
        dat <- ob$corr[[study]]
    }
    `%>%` <- dplyr::`%>%`
    dat <- dat %>%
        dplyr::mutate_at(dplyr::vars(2:ncol(dat)),
                  function(x) x <- ifelse(is.na(x), 0, 1))
    UpSetR::upset(dat, ...)
}

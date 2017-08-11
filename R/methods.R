#' Print method for \link{cmicroRNA} and \link{cTF} objects
#' @param ob A \link{cmicroRNA} or \link{cTF} object such as this returned by
#'    calling \link{cmicroRNA} or \link{cTF}.
#' @param ... Other argument to \code{\link[base]{print}}.
#' @examples
#' # downlaod db file
#' get_db(test = TRUE)
#'
#' # get data for 2 microRNAs in the ACC study
#' dat <- get_mir(c('hsa-let-7b', 'hsa-mir-134'),
#'     study = c('ACC', 'BLCA'),
#'     min_cor = .5)
#'
#' # convert to cmicroRNA object
#' ob <- cmicroRNA(dat)
#' print(ob)
#'
#' @export
print <- function(ob, ...) {
  UseMethod('print')
}

#' @return \code{NULL}
#' @rdname print
#' @method print cmicroRNA
#' @export print cmicroRNA
print.cmicroRNA <- function(ob, ...) {
    p <- paste('A cmicroRNA object: microRNA expression correlations
             in Cancer\n',
               'Contains:\n',
               length(ob$studies), 'Cancer study/ies:', paste(head(ob$studies),
                                                              collapse = ' '),
               '\n',
               length(ob$microRNA), 'microRNA/s:', paste(head(ob$microRNA),
                                                         collapse = ' '),
               '\n',
               length(ob$features), 'features:', paste(head(ob$features),
                                                       collapse = ' '),
               '\n')
    cat(p)
}

#' @return \code{NULL}
#' @rdname print
#' @method print cTF
#' @export print cTF
print.cTF <- function(ob, ...) {
    p <- paste('A cTF object: transcription factor expression correlations
             in Cancer\n',
               'Contains:\n',
               length(ob$studies), 'Cancer study/ies:', paste(head(ob$studies),
                                                              collapse = ' '),
               '\n',
               length(ob$tf), 'Transcription factor/s:', paste(head(ob$TF),
                                                               collapse = ' '),
               '\n',
               length(ob$features), 'features:', paste(head(ob$features),
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
#' @examples
#' # downlaod db file
#' get_db(test = TRUE)
#'
#' # get data for 2 microRNAs in the ACC study
#' dat <- get_mir(c('hsa-let-7b', 'hsa-mir-134'),
#'     study = c('ACC', 'BLCA'),
#'     min_cor = .5)
#'
#' # convert to cmicroRNA object
#' ob <- cmicroRNA(dat)
#' plot(ob, study = 'ACC')
#'
#' @import dplyr tidyr ggplot2
#' @export
plot <- function(ob, study = NULL, ...) {
  UseMethod('plot')
}

#' @return \code{NULL}
#' @rdname plot
#' @method plot cmicroRNA
#' @export plot cmicroRNA
plot.cmicroRNA <- function(ob, study = NULL, ...) {
    if(length(ob$studies) > 1 && length(study) != 1) {
        stop('User shouls provide a singl study to plot.')
    }

    if(is.data.frame(ob$corr)) {
        dat <- ob$corr
    } else {
        dat <- ob$corr[[study]]
    }
    dat <- dat %>%
        gather(mirna_base, cor, -feature) %>%
        na.omit
    gg <- dat %>%
        ggplot(aes(x = mirna_base, y = feature, size = abs(cor))) +
        geom_point()
    return(gg)
}

#' @return \code{NULL}
#' @rdname plot
#' @method plot cTF
#' @export plot cTF
plot.cTF <- function(ob, study = NULL, ...) {
    if(length(ob$studies) > 1 && length(study) != 1) {
        stop('User shouls provide a singl study to plot.')
    }

    if(is.data.frame(ob$corr)) {
        dat <- ob$corr
    } else {
        dat <- ob$corr[[study]]
    }
    dat <- dat %>%
        gather(tf, cor, -feature) %>%
        na.omit
    gg <- dat %>%
        ggplot(aes(x = tf, y = feature, size = abs(cor))) +
        geom_point()
    return(gg)
}

#' Tidy \link{cmicroRNA} objects
#'
#' @inheritParams plot
#'
#' @return A tidy \code{data.frame} of four columns. \code{mirna_base} is the
#' microRNA miRBase IDs, \code{feature} is the features/genes, \code{cor} is
#' the corresponding expression correaltions and \code{study} is TCGA study ID.
#' @examples
#' # downlaod db file
#' get_db(test = TRUE)
#'
#' # get data for 2 microRNAs in the ACC study
#' dat <- get_mir(c('hsa-let-7b', 'hsa-mir-134'),
#'     study = c('ACC', 'BLCA'),
#'     min_cor = .5)
#'
#' # convert to cmicroRNA object
#' ob <- cmicroRNA(dat)
#' tidy(ob) %>%
#'     head
#'
#' @import dplyr tidyr
#' @export
tidy.cmicroRNA <- function(ob) {
    if(length(ob$studies) == 1) {
        dat <- ob$corr %>%
            gather(mirna_base, cor, -feature) %>%
            na.omit %>%
            mutate(study = ob$studies)
    } else {
        dat <- map(ob$corr, function(x) {
            x %>%
                gather(mirna_base, cor, -feature) %>%
                na.omit
        }) %>%
            bind_rows(.id = 'study') %>%
            select(2:3, study)
    }
    return(dat)
}

#' Tidy \link{cTF} objects
#'
#' @inheritParams plot
#'
#' @return A tidy \code{data.frame} of four columns. \code{tf} is the official
#' symbol of the gene containin the transcription factor, \code{feature} is the
#' features/genes, \code{cor} is the corresponding expression correaltions and
#' \code{study} is TCGA study ID.
#' @examples
#' # get db file
#' get_db(test = TRUE)
#'
#' # get data for 2 transcription factors in the ACC study
#' dat <- get_tf(c('AFF4', 'ESR1'),
#'     study = c('ACC', 'BLCA'),
#'     min_cor = .5)
#'
#' # convert to cTF object and plot
#' ob <- cTF(dat)
#' tidy(ob) %>%
#'     head
#'
#' @import dplyr tidyr
#' @export
tidy.cTF <- function(ob) {
    if(length(ob$studies) == 1) {
        dat <- ob$corr %>%
            gather(tf, cor, -feature) %>%
            na.omit %>%
            mutate(study = ob$studies)
    } else {
        dat <- map(ob$corr, function(x) {
            x %>%
                gather(tf, cor, -feature) %>%
                na.omit
        }) %>%
            bind_rows(.id = 'study') %>%
            select(2:3, study)
    }
}

#' Venn Diagram of microRNA correlated features
#'
#' Count and plot the numbers of microRNA correlated features in
#' \code{cmicroRNA} object.
#'
#' @inheritParams plot
#' @return A venn diagram with a ciccle or an ellipses for each microRNA and
#'    the number of correlated features.
#' @examples
#' # get db file
#' get_db(test = TRUE)
#'
#' # get data for 2 microRNAs in the ACC study
#' dat <- get_mir(c('hsa-let-7b', 'hsa-mir-134'),
#'     study = c('ACC', 'BLCA'),
#'     min_cor = .5)
#'
#' # convert to cmicroRNA object and plot
#' ob <- cmicroRNA(dat)
#' venn.diagram(ob, study = 'ACC')
#'
#' @import dplyr tidyr VennDiagram
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
    dat <- dat %>%
        gather(mirna_base, cor, -feature) %>%
        na.omit
    dat <- with(dat, split(feature, mirna_base))
    venn.diagram(dat, ...)
}

#' Venn Diagram of transcription factors correlated features
#'
#' Count and plot the numbers of transcription factors correlated features in
#' \code{cTF} object.
#'
#' @inheritParams plot
#' @return A venn diagram with a ciccle or an ellipses for each transcription
#'    factor and the number of correlated features.
#' @examples
#' # get db file
#' get_db(test = TRUE)
#'
#' # get data for 2 transcription factors in the ACC study
#' dat <- get_tf(c('AFF4', 'ESR1'),
#'     study = c('ACC', 'BLCA'),
#'     min_cor = .5)
#'
#' # convert to cTF object and plot
#' ob <- cTF(dat)
#' venn.diagram(ob, study = 'ACC')
#'
#' @import dplyr tidyr VennDiagram
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
    dat <- dat %>%
        gather(tf, cor, -feature) %>%
        na.omit
    dat <- with(dat, split(feature, tf))
    venn.diagram(dat, ...)
}

#' \code{\link[UpSetR]{upset}} plot of microRNA sets
#'
#' \code{\link[UpSetR]{upset}} of sets of microRNAs and their correlated
#' features in a TCGA study.
#'
#' @inheritParams plot
#' @return An \code{\link[UpSetR]{upset}} plot
#' @examples
#' # get db file
#' get_db(test = TRUE)
#'
#' # get data for 2 microRNAs in the ACC study
#' dat <- get_mir(c('hsa-let-7b', 'hsa-mir-134'),
#'     study = c('ACC', 'BLCA'),
#'     min_cor = .5)
#'
#' # convert to cmicroRNA object and plot
#' ob <- cmicroRNA(dat)
#' upset(ob, study = 'ACC')
#'
#' @import dplyr UpSetR
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

    dat <- dat %>%
        mutate_at(vars(2:ncol(dat)),
                  function(x) x <- ifelse(is.na(x), 0, 1))
    upset(dat, ...)
}

#' \code{\link[UpSetR]{upset}} plot of transcription factors' sets
#'
#' \code{\link[UpSetR]{upset}} of sets of transcription factors and their
#' correlated features in a TCGA study.
#'
#' @inheritParams plot
#' @return An \code{\link[UpSetR]{upset}} plot
#' @examples
#' # get db file
#' get_db(test = TRUE)
#'
#' # get data for 2 transcription factors in the ACC study
#' dat <- get_tf(c('AFF4', 'ESR1'),
#'     study = c('ACC', 'BLCA'),
#'     min_cor = .5)
#'
#' # convert to cTF object and plot
#' ob <- cTF(dat)
#' upset(ob, study = 'ACC')
#'
#' @import dplyr UpSetR
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

    dat <- dat %>%
        mutate_at(vars(2:ncol(dat)),
                  function(x) x <- ifelse(is.na(x), 0, 1))
    upset(dat, ...)
}

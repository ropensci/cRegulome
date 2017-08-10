library(tidyverse)
library(reshape2)
library(UpSetR)
library(VennDiagram)
source('R/get_data.R')
source('R/objects.R')
print.cmicroRNA <- function(ob) {
  p <- paste('A cmicroRNA object: microRNA expression correlations in Cancer\n',
             'Contains:\n',
             length(ob$studies), 'Cancer study/ies:', paste(head(ob$studies), collapse = ' '), '\n',
             length(ob$microRNA), 'microRNA/s:', paste(head(ob$microRNA), collapse = ' '), '\n',
             length(ob$features), 'features:', paste(head(ob$features), collapse = ' '), '\n')
 cat(p)
}

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
}

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

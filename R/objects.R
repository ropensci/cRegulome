cmicroRNA <- function(dat_mir){
  microRNA <- unique(dat_mir$mirna_base)
  features <- unique(dat_mir$feature)
  studies <- unique(dat_mir$study)
  if(length(studies) == 1) {
    corr <- dat_mir %>%
      dcast(feature ~ mirna_base, value.var = 'cor')
  } else {
    corr <- lapply(unique(dat_mir$study),
                   function(x) {
                     dat_mir %>%
                       filter(study == x) %>%
                       dcast(feature ~ mirna_base, value.var = 'cor')
                   })
    names(corr) <- unique(dat_mir$study)
  }
  structure(list(
    microRNA = microRNA,
    features = features,
    studies = studies,
    corr = corr
  ),
  class = 'cmicroRNA')
}

cTF <- function(dat_tf){
  TF <- unique(dat_tf$tf)
  features <- unique(dat_tf$feature)
  studies <- unique(dat_tf$study)
  if(length(studies) == 1) {
    corr <- dat_tf %>%
      dcast(feature ~ tf, value.var = 'cor')
  } else {
    corr <- lapply(unique(dat_tf$study),
                   function(x) {
                     dat_tf %>%
                       filter(study == x) %>%
                       dcast(feature ~ tf, value.var = 'cor')
                   })
    names(corr) <- unique(dat_tf$study)
  }
  structure(list(
    TF = TF,
    features = features,
    studies = studies,
    corr = corr
  ),
  class = 'cTF')
}

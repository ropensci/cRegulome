get_db <- function(...) {
  # db file url
  url <- 'https://www.dropbox.com/s/unaj94rnk0n4fl5/cRegulome.db.gz?raw=1'

  # check url exists
  if(!RCurl::url.exists(url)) {
    stop("URL doesn't exist.")
  }

  # download file
  if(file.exists('cRegulome.db.gz')) {
    message('File already exists in the current directory.')
  } else {
    tryCatch(download.file(url, destfile = 'cRegulome.db.gz'),
             error = function() {
               message('File download failed.')
               return(NA)
             })
  }
}

get_mir <- function(mir, study = NULL, min_cor = NULL, max_num = NULL, targets_only = FALSE) {
  # check db file exist in the current directory
  if(!file.exists('cRegulome.db')){
    stop("Database file doesn't exist in the current directory.")
  }

  # connect to db
  db <- DBI::dbConnect(RSQLite::SQLite(),
                       'cRegulome.db')

  # unpack filters and check types
  table <- 'cor_mir'

  if(is.null(mir)) {
    stop("User should provide at least one microRNA ID")
  } else if(!is.character(mir)) {
    stop("mir should be a character vector")
  } else {
    mir <- as.character(mir)
  }

  if(is.null(study)) {
    study <- DBI::dbListFields(db, table)[-1:-2]
  } else if(!is.character(study)){
    stop("Study should be a character vector")
  } else {
    study <- as.character(study)
  }

  if(is.null(min_cor)) {
    min_cor <- 0
  } else if(!is.numeric(min_cor) || min_cor > 1 || min_cor < 0) {
    stop("min_cor should be a numeric between 0 and 1.")
  }
  else {
    min_cor <- as.numeric(min_cor)
  }

  # get main data by applying filters and tidy
  dat <- db %>%
    tbl(table) %>%
    select(mirna_base, feature, study) %>%
    filter(mirna_base %in% mir) %>%
    collect %>%
    tidyr::gather(study, cor, -mirna_base, -feature) %>%
    mutate(cor = cor/100) %>%
    filter(abs(cor) > min_cor) %>%
    arrange(desc(abs(cor))) %>%
    na.omit

  # apply targets only filters when TRUE
  if(targets_only == TRUE) {
    # subset targets
    targets <- db %>%
      tbl('targets_mir') %>%
      filter(mirna_base %in% mir) %>%
      collect

    # subset main data to targets only
    dat <- right_join(dat, targets) %>%
      na.omit
  }

  # subset to max_num
  if(is.null(max_num)) {
  } else if(is.integer(max_num) || max_num < 0) {
    stop("max_num should be an integer between 1 and Inf.")
  } else {
    dat <- dat %>%
      group_by(mirna_base, study) %>%
      slice(1:max_num)
  }

  # disconnect and return dat
  DBI::dbDisconnect(db)
  return(dat)
}

get_tf <- function(tf, study = NULL, min_cor = NULL, max_num = NULL, targets_only = FALSE) {
  # check db file exist in the current directory
  if(!file.exists('cRegulome.db')){
    stop("Database file doesn't exist in the current directory.")
  }

  # connect to db
  db <- DBI::dbConnect(RSQLite::SQLite(),
                       'cRegulome.db')

  # unpack filters and check types
  table <- 'cor_tf'

  if(is.null(tf)) {
    stop("User should provide at least one microRNA ID")
  } else if(!is.character(tf)) {
    stop("tf should be a character vector")
  } else {
    tf_id <- as.character(tf)
  }

  if(is.null(study)) {
    studies <- DBI::dbListFields(db, table)[-1:-2]
  } else if(!is.character(study)){
    stop("Study should be a character vector")
  } else {
    study <- as.character(study)
  }

  if(is.null(min_cor)) {
    min_cor <- 0
  } else if(!is.numeric(min_cor) || min_cor > 1 || min_cor < 0) {
    stop("min_cor should be a numeric between 0 and 1.")
  }
  else {
    min_cor <- as.numeric(min_cor)
  }

  # get main data by applying filters and tidy
  dat <- db %>%
    tbl(table) %>%
    select(tf, feature, study) %>%
    filter(tf %in% tf_id) %>%
    collect %>%
    tidyr::gather(study, cor, -tf, -feature) %>%
    mutate(cor = cor/100) %>%
    filter(abs(cor) > min_cor) %>%
    arrange(desc(abs(cor))) %>%
    na.omit

  # apply targets only filters when TRUE
  if(targets_only == TRUE) {
    # subset targets
    targets <- db %>%
      tbl('targets_tf') %>%
      filter(tf %in% tf_id) %>%
      collect

    # subset main data to targets only
    dat <- right_join(dat, targets) %>%
      na.omit
  }

  # subset to max_num
  if(is.null(max_num)) {
  } else if(is.integer(max_num) || max_num < 0) {
    stop("max_num should be an integer between 1 and Inf.")
  } else {
    dat <- dat %>%
      group_by(tf, study) %>%
      slice(1:max_num)
  }

  # disconnect and return dat
  DBI::dbDisconnect(db)
  return(dat)
}



context('Downlaoding database and extracting data')
test_that('get_db downlaods db file succussfully', {
  library(RCurl)
  get_db(test = TRUE)
  expect_true(file.exists('cRegulome.db.gz'))
  if(file.exists('cRegulome.db.gz')) unlink('cRegulome.db.gz')
})

test_that('get_mir extracts data properly', {
  # load libraries
  library(RCurl)
  library(R.utils)

  # get error when file doesn't exist
  expect_error(get_mir('hsa-let-7b', 'ACC'))

  # get db file
  get_db(test = TRUE)
  gunzip('cRegulome.db.gz')

  # falty arguments
  expect_error(get_mir(study = 'ACC'))
  expect_error(get_mir('hsa-let-7b', min_cor = -1))

  # chech output is a data.frame and nrows > 1
  dat <- get_mir('hsa-let-7b', 'ACC')
  expect_true(is.data.frame(dat))
  expect_true(nrow(dat) > 1)

  # when targets_only
  dat_targets <- get_mir('hsa-let-7b', 'ACC', targets_only = TRUE)
  expect_true(nrow(dat) > nrow(dat_targets))

  # clean up
  if(file.exists('cRegulome.db')) unlink('cRegulome.db')
})

test_that('objects cmicroRNA are constructed properly', {
  # load libraries
  library(RCurl)
  library(R.utils)

  # get error when file doesn't exist
  expect_error(get_mir('hsa-let-7b', 'ACC'))

  # get db file
  get_db(test = TRUE)
  gunzip('cRegulome.db.gz')
  dat <- get_mir('hsa-let-7b', 'ACC', min_cor = .3, max_num = 5)
  cmir <- cmicroRNA(dat)

  # check class and names
  expect_identical(class(cmir), 'cmicroRNA')
  expect_identical(names(cmir), c("microRNA", "features", "studies", "corr"))

  # check components
  expect_identical(cmir$microRNA, 'hsa-let-7b')
  expect_identical(cmir$studies, 'ACC')
  expect_equal(length(cmir$features), 5)
  expect_true(min(abs(cmir$corr$`hsa-let-7b`)) > .3)

  # clean up
  if(file.exists('cRegulome.db')) unlink('cRegulome.db')
})

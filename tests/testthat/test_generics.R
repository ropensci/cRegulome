context('Downlaoding database and extracting data')


test_that('get_db downlaods db file succussfully', {
    get_db(test = TRUE)
    expect_true(file.exists('cRegulome.db.gz'))
    })

# decompress database file
R.utils::gunzip('cRegulome.db.gz')

test_that('get_db stop when db file exist in the current directory', {
    expect_message(get_db(test = TRUE))
})

# connect to the db file
conn <- DBI::dbConnect(RSQLite::SQLite(), 'cRegulome.db')

test_that('get_mir Faulty arguments', {
    expect_error(get_mir(conn, mir = NULL))
    expect_error(get_mir(conn, mir = list('hsa-let-7b')))
    expect_error(get_mir(conn, study = 'ACC'))
    expect_error(get_mir(conn, mir = 'hsa-let-7b', study = 11))
    expect_error(get_mir(conn, mir = 'hsa-let-7b', min_cor = -1))
    expect_error(get_mir(conn, mir = 'hsa-let-7b', max_num = -1))
    })

# get_mir with and without targets_only
dat <- get_mir(conn,
               mir = 'hsa-let-7b',
               study = 'ACC')
dat_targets <- get_mir(conn,
                       mir = 'hsa-let-7b',
                       study = 'ACC',
                       targets_only = TRUE)

test_that('get_mir extracts data properly', {
    expect_true(is.data.frame(dat))
    expect_true(nrow(dat) > 1)
    expect_true(nrow(dat) > nrow(dat_targets))
    })

# make a cmircroRNA object
dat <- get_mir(conn,
               mir = 'hsa-let-7b',
               study = 'ACC',
               min_cor = .3,
               max_num = 5)
cmir <- cmicroRNA(dat)

test_that('objects cmicroRNA are constructed properly', {
    expect_s3_class(cmir, 'cmicroRNA')
    expect_identical(names(cmir), c("microRNA", "features", "studies", "corr"))
    expect_identical(cmir$microRNA, 'hsa-let-7b')
    expect_identical(cmir$studies, 'ACC')
    expect_equal(length(cmir$features), 5)
    expect_true(min(abs(cmir$corr$`hsa-let-7b`)) > .3)
    })

dat <- get_mir(conn,
               mir = 'hsa-let-7b',
               study = c('ACC', 'BLCA'),
               min_cor = .3)
cmir <- cmicroRNA(dat)

test_that('object cmicroRNA with multiple studies', {
  expect_s3_class(cmir, 'cmicroRNA')
  expect_identical(names(cmir), c("microRNA", "features", "studies", "corr"))
  expect_identical(cmir$studies, c('ACC', 'BLCA'))
})

test_that('get_tf Faulty arguments', {
  expect_error(get_tf(conn, tf = NULL))
  expect_error(get_tf(conn, tf = list('AFF4')))
  expect_error(get_tf(conn, study = 'ACC'))
  expect_error(get_tf(conn, tf = 'AFF4', study = 11))
  expect_error(get_tf(conn, tf = 'AFF4', min_cor = -1))
  expect_error(get_tf(conn, tf = 'AFF4', max_num = -1))
})

# get_tf with and without targets_only
dat <- get_tf(conn,
               tf = 'AFF4',
               study = 'ACC')
dat_targets <- get_tf(conn,
                       tf = 'AFF',
                       study = 'ACC',
                       targets_only = TRUE)

test_that('get_tf extracts data properly', {
  expect_true(is.data.frame(dat))
  expect_true(nrow(dat) > 1)
  expect_true(nrow(dat) > nrow(dat_targets))
})

# make a cmircroRNA object
dat <- get_tf(conn,
               tf = 'AFF4',
               study = 'ACC',
               min_cor = .3,
               max_num = 5)
ctf <- cTF(dat)

test_that('objects cTF are constructed properly', {
  expect_s3_class(ctf, 'cTF')
  expect_identical(names(ctf), c("TF", "features", "studies", "corr"))
  expect_identical(ctf$TF, 'AFF4')
  expect_identical(ctf$studies, 'ACC')
  expect_equal(length(ctf$features), 5)
  expect_true(min(abs(ctf$corr$`AFF4`)) > .3)
})

dat <- get_tf(conn,
              tf = c('AFF4', 'AR'),
              study = c('ACC', 'BLCA'),
              min_cor = .3)
ctf <- cTF(dat)

test_that('object cTF with multiple studies', {
  expect_s3_class(ctf, 'cTF')
  expect_identical(names(ctf), c("TF", "features", "studies", "corr"))
  expect_identical(ctf$studies, c('ACC', 'BLCA'))
})

# clean up
DBI::dbDisconnect(conn)
if(file.exists('cRegulome.db')) unlink('cRegulome.db')

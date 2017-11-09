context('Downlaoding database and extracting data')

# if(file.exists('cRegulome.db')) unlink('cRegulome.db')

# test_that('get_db downlaods db file succussfully', {
#     get_db(test = TRUE)
#     expect_true(file.exists('cRegulome.db.gz'))
# })

# # decompress database file
# R.utils::gunzip('cRegulome.db.gz')

# test_that('get_db stop when db file exist in the current directory', {
#    expect_message(get_db(test = TRUE))
# })
# if(file.exists('cRegulome.db')) unlink('cRegulome.db')

# connect to the db file
fl <- system.file('extdata', 'cRegulome.db', package = 'cRegulome')
conn <- dbConnect(SQLite(), fl)

test_that('get_mir Faulty arguments', {
    expect_error(get_mir(conn, mir = NULL),
                 'User should provide at least one microRNA ID')
    expect_error(get_mir(conn, mir = list('hsa-let-7g')),
                 'mir should be a character vector')
    expect_error(get_mir(conn, study = 'STES'),
                 'argument "mir" is missing, with no default')
    expect_error(get_mir(conn, mir = 'hsa-let-7g', study = 11),
                 'Study should be a character vector')
    expect_error(get_mir(conn, mir = 'hsa-let-7g', min_abs_cor = -1),
                 'min_abs_cor should be a numeric between 0 and 1.')
    expect_error(get_mir(conn, mir = 'hsa-let-7g', max_num = -1),
                 'max_num should be a positive integer.')
    expect_error(get_mir(conn, mir = 'hsa-let-7g', max_num = 0),
                 'max_num should be a positive integer.')
    expect_error(get_mir(conn, mir = 'hsa-let-7g', max_num = .5),
                 'max_num should be a positive integer.')
    })

test_that('get_mir extracts data properly', {
    dat <- get_mir(conn,
                   mir = 'hsa-let-7g',
                   study = 'STES')
    
    expect_true(is.data.frame(dat))
    expect_true(nrow(dat) > 1)
    })

test_that('objects cmicroRNA are constructed properly', {
    dat <- get_mir(conn,
                   mir = 'hsa-let-7g',
                   study = 'STES',
                   min_abs_cor = .3,
                   max_num = 5)
    cmir <- cmicroRNA(dat)
    
    expect_s3_class(cmir, 'cmicroRNA')
    expect_identical(names(cmir), c("microRNA", "features", "studies", "corr"))
    expect_identical(cmir$microRNA, 'hsa-let-7g')
    expect_identical(cmir$studies, 'STES')
    expect_equal(length(cmir$features), 5)
    expect_true(min(abs(cmir$corr$`hsa-let-7g`)) > .3)
    })

test_that('get_tf Faulty arguments', {
  expect_error(get_tf(conn, tf = NULL), 
               'User should provide at least one TF ID')
  expect_error(get_tf(conn, tf = list('LEF1')),
               'tf should be a character vector')
  expect_error(get_tf(conn, study = 'STES'),
               'argument "tf" is missing, with no default')
  expect_error(get_tf(conn, tf = 'LEF1', study = 11),
               'Study should be a character vector')
  expect_error(get_tf(conn, tf = 'LEF1', min_abs_cor = -1),
               'min_abs_cor should be a numeric between 0 and 1.')
  expect_error(get_tf(conn, tf = 'LEF1', max_num = -1),
               'max_num should be a positive integer.')
  expect_error(get_tf(conn, tf = 'LEF1', max_num = 0),
               'max_num should be a positive integer.')
  expect_error(get_tf(conn, tf = 'LEF1', max_num = .5),
               'max_num should be a positive integer.')
})

test_that('get_tf extracts data properly', {
    dat <- get_tf(conn,
                  tf = 'LEF1',
                  study = 'STES*')
    
  expect_true(is.data.frame(dat))
  expect_true(nrow(dat) > 1)
})

# make a cmircroRNA object


test_that('objects cTF are constructed properly', {
    dat <- get_tf(conn,
                  tf = 'LEF1',
                  study = 'STES*',
                  min_abs_cor = .3,
                  max_num = 5)
    ctf <- cTF(dat)
    
  expect_s3_class(ctf, 'cTF')
  expect_identical(names(ctf), c("TF", "features", "studies", "corr"))
  expect_identical(ctf$TF, 'LEF1')
  expect_identical(ctf$studies, 'STES*')
  expect_equal(length(ctf$features), 5)
  expect_true(min(abs(ctf$corr$`LEF1`)) > .3)
})

# clean up
dbDisconnect(conn)
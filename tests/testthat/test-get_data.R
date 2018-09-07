context('get_data')

# remove db file if exists
if(file.exists('cRegulome.db')) {
    unlink('cRegulome.db')
}

test_that('get_db downlaods db file succussfully', {
    get_db(test = TRUE)
    fl <- paste(tempdir(), 'cRegulome.db', sep = '/')
    
    expect_true(file.exists(fl))
    
    if(file.exists(fl)) unlink(fl)
})

test_that('get_db downloads db file to certain location', {
    get_db(test = TRUE,
           destfile = './cRegulome.db')
    
    fl <- './cRegulome.db'
    expect_true(file.exists(fl))
    
    if(file.exists(fl)) unlink(fl)
})

# connect to the db file
fl <- system.file('extdata', 'cRegulome.db', package = 'cRegulome')
conn <- RSQLite::dbConnect(RSQLite::SQLite(), fl)

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

test_that('get_tf extracts data properly', {
    dat <- get_tf(conn,
                  tf = 'LEF1',
                  study = 'STES*')
    
    expect_true(is.data.frame(dat))
    expect_true(nrow(dat) > 1)
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

# clean up
RSQLite::dbDisconnect(conn)
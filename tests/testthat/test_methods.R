context('Testing methods output')

# connect to the db file
fl <- system.file('extdata', 'cRegulome.db', package = 'cRegulome')
conn <- dbConnect(SQLite(), fl)

test_that('cor_tidy output a data.frame identical to output of get_mir', {
    dat <- get_mir(conn,
                   mir = 'hsa-let-7g',
                   study = 'STES',
                   min_abs_cor = .3,
                   max_num = 5)
    
    cmir <- cmicroRNA(dat)
    tidy_cmir <- cor_tidy(cmir)
    
    expect_equal(dat, tidy_cmir)
})


test_that('tidy output a data.frame identical to output of get_tf', {
    dat <- get_tf(conn,
                  tf = 'LEF1',
                  study = 'STES*',
                  min_abs_cor = .3,
                  max_num = 5)
    
    ctf <- cTF(dat)
    tidy_tf <- cor_tidy(ctf)
    
    expect_equal(dat, tidy_tf)
})


test_that("cor_igraph retruns a the proper object", {
    dat <- get_mir(conn,
                   mir = 'hsa-let-7g',
                   study = 'STES',
                   min_abs_cor = .3,
                   max_num = 5)
    
    cmir <- cmicroRNA(dat)
    
    g <- cor_igraph(cmir)
    expect_s3_class(g, 'igraph')
})

# clean up
dbDisconnect(conn)

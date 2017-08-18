context('Downlaoding database and extracting data')

test_that('get_db downlaods db file succussfully', {
    get_db(test = TRUE)
    expect_true(file.exists('cRegulome.db.gz'))
    })

# decompress database file
R.utils::gunzip('cRegulome.db.gz')
conn <- DBI::dbConnect(RSQLite::SQLite(), 'cRegulome.db')

test_that('Faulty arguments', {
    expect_error(get_mir(conn,
                         study = 'ACC'))
    expect_error(get_mir(conn,
                         mir = 'hsa-let-7b',
                         min_cor = -1))
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
    expect_identical(class(cmir), 'cmicroRNA')
    expect_identical(names(cmir), c("microRNA", "features", "studies", "corr"))
    expect_identical(cmir$microRNA, 'hsa-let-7b')
    expect_identical(cmir$studies, 'ACC')
    expect_equal(length(cmir$features), 5)
    expect_true(min(abs(cmir$corr$`hsa-let-7b`)) > .3)
    })

# clean up
DBI::dbDisconnect(conn)
if(file.exists('cRegulome.db')) unlink('cRegulome.db')

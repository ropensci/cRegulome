context('Testing methods output')

# get db file and make object
get_db(test = TRUE)
R.utils::gunzip('cRegulome.db.gz')

conn <- DBI::dbConnect(RSQLite::SQLite(), 'cRegulome.db')
dat <- get_mir(conn,
               mir = 'hsa-let-7b',
               study = 'ACC',
               min_cor = .3,
               max_num = 5)

cmir <- cmicroRNA(dat)
tidy_cmir <- tidy(cmir)

test_that('tidy output a data.frame identical to output of get_mir', {
    expect_equal(dat, tidy_cmir)
    })

dat <- get_tf(conn,
               tf = 'AFF4',
               study = 'ACC',
               min_cor = .3,
               max_num = 5)

ctf <- cTF(dat)
tidy_tf <- tidy(ctf)

test_that('tidy output a data.frame identical to output of get_tf', {
  expect_equal(dat, tidy_tf)
})

# clean up
DBI::dbDisconnect(conn)
if(file.exists('cRegulome.db')) unlink('cRegulome.db')

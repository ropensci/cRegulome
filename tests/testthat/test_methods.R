context('Testing methods output')

# get db file and make object
get_db(test = TRUE)
R.utils::gunzip('cRegulome.db.gz')

conn <- DBI::dbConnect(RSQLite::SQLite(), 'cRegulome.db')
dat1 <- get_mir(conn,
               mir = 'hsa-let-7b',
               study = 'ACC',
               min_abs_cor = .3,
               max_num = 5)

cmir1 <- cmicroRNA(dat1)
tidy_cmir1 <- tidy(cmir1)

dat2 <- get_mir(conn,
                mir = 'hsa-let-7b',
                study = c('ACC', 'BLCA'),
                min_abs_cor = .3,
                max_num = 5)

cmir2 <- cmicroRNA(dat2)
tidy_cmir2 <- tidy(cmir2)

test_that('tidy output a data.frame identical to output of get_mir', {
  expect_equal(dat1, tidy_cmir1)
  expect_equal(dat2, tidy_cmir2)
})


dat1 <- get_tf(conn,
               tf = 'AFF4',
               study = 'ACC',
               min_abs_cor = .3,
               max_num = 5)

ctf1 <- cTF(dat1)
tidy_tf1 <- tidy(ctf1)

dat2 <- get_tf(conn,
               tf = 'AFF4',
               study = c('ACC', 'BLCA'),
               min_abs_cor = .3,
               max_num = 5)

ctf2 <- cTF(dat2)
tidy_tf2 <- tidy(ctf2)

test_that('tidy output a data.frame identical to output of get_tf', {
  expect_equal(dat1, tidy_tf1)
  expect_equal(dat2, tidy_tf2)
})

# clean up
DBI::dbDisconnect(conn)
if(file.exists('cRegulome.db')) unlink('cRegulome.db')

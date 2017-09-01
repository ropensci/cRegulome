# decompress database file
get_db(test = TRUE)
R.utils::gunzip('cRegulome.db.gz')
conn <- DBI::dbConnect(RSQLite::SQLite(), 'cRegulome.db')

conn <- DBI::dbConnect(RSQLite::SQLite(), 'cRegulome.db')
dat <- get_mir(conn,
               mir = c('hsa-let-7b', 'hsa-mir-134'),
               study = 'ACC',
               min_cor = .3)

cmir <- cmicroRNA(dat)

test_that('plot cmicroRNA objects', {
  gg <- plot(cmir)  
  expect_identical(class(gg), c('gg', 'ggplot'))
  expect_equal(dim(gg$data), c(nrow(dat), 5))
  
  expect_identical(class(gg$layers[[1]]$geom),
                   c("GeomPoint", "Geom", "ggproto"))
  expect_equal(gg$labels, list(size = 'Correlation',
                               colour = 'Direction',
                               x = 'mirna_base',
                               y = 'feature'))
})

test_that('upset cmicroRNA objects', {
  gg <- joy(cmir)  
  expect_identical(class(gg), c('gg', 'ggplot'))
  expect_equal(dim(gg$data), c(nrow(dat), 3))
  
  expect_identical(class(gg$layers[[1]]$geom),
                   c("GeomJoy", "GeomRidgeline", "Geom", "ggproto"))
  expect_equal(gg$labels, list(x = 'cor',
                               y = 'mirna_base',
                               height = 'density'))
})


dat <- get_tf(conn,
               tf = c('AFF4', 'AR'),
               study = 'ACC',
               min_cor = .3)

ctf <- cTF(dat)

test_that('plot cTF objects', {
  gg <- plot(ctf)  
  expect_identical(class(gg), c('gg', 'ggplot'))
  expect_equal(dim(gg$data), c(nrow(dat), 5))
  
  expect_identical(class(gg$layers[[1]]$geom),
                   c("GeomPoint", "Geom", "ggproto"))
  expect_equal(gg$labels, list(size = 'Correlation',
                               colour = 'Direction',
                               x = 'tf',
                               y = 'feature'))
})

test_that('upset cTF objects', {
  gg <- joy(ctf)  
  expect_identical(class(gg), c('gg', 'ggplot'))
  expect_equal(dim(gg$data), c(nrow(dat), 3))
  
  expect_identical(class(gg$layers[[1]]$geom),
                   c("GeomJoy", "GeomRidgeline", "Geom", "ggproto"))
  expect_equal(gg$labels, list(x = 'cor',
                               y = 'tf',
                               height = 'density'))
})

# clean up
DBI::dbDisconnect(conn)
if(file.exists('cRegulome.db')) unlink('cRegulome.db')


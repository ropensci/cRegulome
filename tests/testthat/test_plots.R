# decompress database file
get_db(test = TRUE)
R.utils::gunzip('cRegulome.db.gz')
conn <- DBI::dbConnect(RSQLite::SQLite(), 'cRegulome.db')
dat1 <- get_mir(conn,
                mir = c('hsa-let-7b', 'hsa-mir-134'),
                study = 'ACC',
                min_abs_cor = .3)
dat2 <- get_mir(conn,
               mir = c('hsa-let-7b', 'hsa-mir-134'),
               study = c('ACC', 'BLCA'),
               min_abs_cor = .3)

cmir1 <- cmicroRNA(dat1)
cmir2 <- cmicroRNA(dat2)

DBI::dbDisconnect(conn)

test_that('cor_plot cmicroRNA objects', {
    gg <- cor_plot(cmir1)  
    expect_s3_class(gg, c('gg', 'ggplot'))
    expect_equal(dim(gg$data), c(nrow(dat1), 5))
    
    expect_s3_class(gg$layers[[1]]$geom,
                     c("GeomPoint", "Geom", "ggproto"))
    expect_equal(gg$labels, list(size = 'Correlation',
                                 colour = 'Direction',
                                 x = 'mirna_base',
                                 y = 'feature'))
    expect_error(cor_plot(cmir2))
    gg <- cor_plot(cmir2, study = 'ACC')
    expect_s3_class(gg, c('gg', 'ggplot'))
})

test_that('cor_joy cmicroRNA objects', {
    gg <- cor_joy(cmir1)  
    expect_s3_class(gg, c('gg', 'ggplot'))
    expect_equal(dim(gg$data), c(nrow(dat1), 3))
    
    expect_s3_class(gg$layers[[1]]$geom,
                     c("GeomDensityRidges", "GeomRidgeline", "Geom", "ggproto"))
    expect_equal(gg$labels, list(x = 'cor',
                                 y = 'mirna_base',
                                 height = 'density'))
    expect_error(cor_joy(cmir2))
    gg <- cor_joy(cmir2, study = 'ACC')
    expect_s3_class(gg, c('gg', 'ggplot'))
})

conn <- DBI::dbConnect(RSQLite::SQLite(), 'cRegulome.db')
dat1 <- get_tf(conn,
               tf = c('AFF4', 'AR'),
               study = 'ACC',
               min_abs_cor = .3)
dat2 <- get_tf(conn,
               tf = c('AFF4', 'AR'),
               study = c('ACC', 'BLCA'),
               min_abs_cor = .3)

ctf1 <- cTF(dat1)
ctf2 <- cTF(dat2)

DBI::dbDisconnect(conn)

test_that('cor_plot cTF objects', {
    gg <- cor_plot(ctf1)  
    expect_s3_class(gg, c('gg', 'ggplot'))
    expect_equal(dim(gg$data), c(nrow(dat1), 5))
    
    expect_s3_class(gg$layers[[1]]$geom,
                     c("GeomPoint", "Geom", "ggproto"))
    expect_equal(gg$labels, list(size = 'Correlation',
                                 colour = 'Direction',
                                 x = 'tf',
                                 y = 'feature'))
    expect_error(cor_plot(ctf2))
    gg <- cor_plot(ctf2, study = 'ACC')
    expect_s3_class(gg, c('gg', 'ggplot'))
})

test_that('cor_joy cTF objects', {
    gg <- cor_joy(ctf1)  
    expect_s3_class(gg, c('gg', 'ggplot'))
    expect_equal(dim(gg$data), c(nrow(dat1), 3))
    
    expect_s3_class(gg$layers[[1]]$geom,
                     c("GeomDensityRidges", "GeomRidgeline", "Geom", "ggproto"))
    expect_equal(gg$labels, list(x = 'cor',
                                 y = 'tf',
                                 height = 'density'))
    expect_error(cor_joy(ctf2))
    gg <- cor_joy(ctf2, study = 'ACC')
    expect_s3_class(gg, c('gg', 'ggplot'))
})

# clean up
if(file.exists('cRegulome.db')) unlink('cRegulome.db')


context("Testing plots")

# connect to the db file
fl <- system.file('extdata', 'cRegulome.db', package = 'cRegulome')
conn <- dbConnect(SQLite(), fl)
dat <- get_mir(conn,
                mir = 'hsa-let-7g',
                study = 'STES',
                min_abs_cor = .3)
cmir <- cmicroRNA(dat)

test_that('cor_plot cmicroRNA objects', {
    gg <- cor_plot(cmir)  
    expect_s3_class(gg, c('gg', 'ggplot'))
    expect_equal(dim(gg$data), c(nrow(dat), 5))
    
    expect_s3_class(gg$layers[[1]]$geom,
                     c("GeomPoint", "Geom", "ggproto"))
    expect_equal(gg$labels, list(size = 'Correlation',
                                 colour = 'Direction',
                                 x = 'mirna_base',
                                 y = 'feature'))
})

test_that('cor_joy cmicroRNA objects', {
    gg <- cor_joy(cmir)  
    expect_s3_class(gg, c('gg', 'ggplot'))
    expect_equal(dim(gg$data), c(nrow(dat), 3))
    
    expect_s3_class(gg$layers[[1]]$geom,
                     c("GeomDensityRidges", "GeomRidgeline", "Geom", "ggproto"))
    expect_equal(gg$labels, list(x = 'cor',
                                 y = 'mirna_base',
                                 height = 'density'))
})

dat <- get_tf(conn,
               tf = 'LEF1',
               study = 'STES*',
               min_abs_cor = .3)

ctf <- cTF(dat)


test_that('cor_plot cTF objects', {
    gg <- cor_plot(ctf)  
    expect_s3_class(gg, c('gg', 'ggplot'))
    expect_equal(dim(gg$data), c(nrow(dat), 5))
    
    expect_s3_class(gg$layers[[1]]$geom,
                     c("GeomPoint", "Geom", "ggproto"))
    expect_equal(gg$labels, list(size = 'Correlation',
                                 colour = 'Direction',
                                 x = 'tf',
                                 y = 'feature'))
})

test_that('cor_joy cTF objects', {
    gg <- cor_joy(ctf)  
    expect_s3_class(gg, c('gg', 'ggplot'))
    expect_equal(dim(gg$data), c(nrow(dat), 3))
    
    expect_s3_class(gg$layers[[1]]$geom,
                     c("GeomDensityRidges", "GeomRidgeline", "Geom", "ggproto"))
    expect_equal(gg$labels, list(x = 'cor',
                                 y = 'tf',
                                 height = 'density'))
})

# clean up
dbDisconnect(conn)


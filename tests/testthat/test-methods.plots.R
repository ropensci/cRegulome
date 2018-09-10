context("methods.plots")

# locate the testset file and connect
fl <- system.file('extdata', 'cRegulome.db', package = 'cRegulome')
conn <- RSQLite::dbConnect(RSQLite::SQLite(), fl)

# enter a custom query with different arguments
dat <- get_mir(conn,
               mir = c('hsa-let-7g', 'hsa-let-7i'),
               max_num = 5,
               study = 'STES')

# make a cmicroRNA object   
cmir <- cmicroRNA(dat)



# clean up
RSQLite::dbDisconnect(conn)

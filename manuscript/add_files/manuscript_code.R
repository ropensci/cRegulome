library(cRegulome)
library(R.utils)
library(DBI)
library(RSQLite)
library(ggplot2)

# only in final script
# get_db()
# gunzip('cRegulome.db.gz')
# conn <- dbConnect(SQLite(), 'cRegulome.db')

db_file <- system.file("extdata", "cRegulome.db", package = "cRegulome")
conn <- DBI::dbConnect(RSQLite::SQLite(), db_file)

dat <- get_mir(conn,
               mir = c('hsa-let-7b', 'hsa-mir-134'),
               study = c('ACC', 'BLCA'),
               targets_only = TRUE)

dbDisconnect(conn)

ob <- cmicroRNA(dat)
dat <- tidy(ob)
dat <- subset(dat, abs(dat$cor) > .3)
sub_ob <- cmicroRNA(dat)

print(sub_ob)

plot(sub_ob, study = 'ACC') +
  labs(x = '', y = '') +
  theme(legend.position = 'top',
        legend.box = 'vertical')
ggsave(filename = '~/microRNA/cRegulome_manuscript/figures/fig1a.png',
       dpi = 600, width = 4, height = 8, units = 'in')
  
png(filename = '~/microRNA/cRegulome_manuscript/figures/fig1b.png',
    res = 300, width = 4, height = 4, units = 'in', pointsize = 10)
hist(ob, study = 'ACC', breaks = 100, main = '', xlab = 'Correlation')
dev.off()

joy(ob, study = 'ACC') +
  labs(x = 'Correlation', y = '')
ggsave(filename = '~/microRNA/cRegulome_manuscript/figures/fig1c.png',
       dpi = 600, width = 4, height = 4, units = 'in')  

png(filename = '~/microRNA/cRegulome_manuscript/figures/fig2a.png',
    res = 300, width = 4, height = 4, units = 'in', pointsize = 10)
venn.diagram(ob, study = 'ACC', cat.default.pos = 'text')
dev.off()

png(filename = '~/microRNA/cRegulome_manuscript/figures/fig2b.png',
    res = 300, width = 4, height = 4, units = 'in', pointsize = 10)
upset(ob, study = 'ACC')
dev.off()

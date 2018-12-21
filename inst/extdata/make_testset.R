require(RSQLite)
require(dplyr)
library(cRegulome)
library(reshape2)

db <- dbConnect(SQLite(), '../cregart/cRegulome.db')
mirna <- c("hsa-miR-18b*", "hsa-miR-409-3p", "hsa-let-7g", "hsa-miR-30c", "hsa-miR-30a", "hsa-miR-150", "hsa-miR-1237",
           "hsa-miR-200a", "hsa-let-7f", "hsa-miR-92b", "hsa-miR-497", "hsa-miR-99a", "hsa-miR-1247", "hsa-miR-200b*",
           "hsa-miR-940", "hsa-miR-99b", "hsa-miR-30b", "hsa-miR-1227", "hsa-miR-501-3p", "hsa-miR-324-3p", "hsa-miR-574-5p",
           "hsa-miR-486-5p", "hsa-miR-151-5p", "hsa-miR-29b-2*", "hsa-miR-532-5p", "hsa-miR-342-3p", "hsa-miR-339-3p", "hsa-miR-31",
           "hsa-miR-27b", "hsa-miR-877", "hsa-miR-146a", "hsa-miR-126", "hsa-miR-28-5p", "hsa-miR-125a-5p", "hsa-miR-425*", "hsa-miR-181b", 
           "hsa-miR-193b", "hsa-miR-1225-3p", "hsa-miR-130b", "hsa-miR-652", "hsa-miR-375", "hsa-miR-487b", "hsa-miR-361-5p", "hsa-miR-378",
           "hsa-let-7i", "hsa-miR-34a", "hsa-miR-625", "hsa-miR-30d", "hsa-miR-362-5p", "hsa-miR-92a", "hsa-miR-149", "hsa-miR-423-3p", 
           "hsa-miR-222", "hsa-miR-665", "hsa-miR-221", "hsa-miR-140-3p", "hsa-miR-424*", "hsa-miR-432", "hsa-miR-200b", "hsa-miR-151-3p",
           "hsa-miR-192", "hsa-miR-324-5p", "hsa-miR-195", "hsa-miR-423-5p", "hsa-miR-378*", "hsa-miR-27a", "hsa-miR-502-5p", "hsa-miR-585",
           "hsa-miR-488*", "hsa-miR-424", "hsa-miR-1274b", "hsa-miR-216b", "hsa-miR-101", "hsa-miR-362-3p", "hsa-miR-509-5p", "hsa-miR-511", 
           "hsa-miR-182*", "hsa-miR-1253", "hsa-miR-542-3p", "hsa-miR-1260", "hsa-miR-331-3p", "hsa-miR-920", "hsa-miR-217", "hsa-miR-923", 
           "hsa-miR-155*", "hsa-miR-33a", "hsa-miR-888", "hsa-miR-566", "hsa-miR-431", "hsa-miR-518a-3p", "hsa-miR-1826", "hsa-miR-202", 
           "hsa-miR-147b")

cor_mir <- get_mir(db,
                   mir = tolower(mirna),
                   study = 'STES',
                   targets_only = TRUE)
cor_mir <- mutate(cor_mir, cor = cor * 100) %>%
    dcast(mirna_base + feature ~ study, value.var = 'cor')
targets_mir <- select(cor_mir, mirna_base, feature) %>% unique()


tf_id <- c("ETV4", "LEF1", "MYB", "MYBL2", "TFAP2A")

cor_tf <- get_tf(db,
                 tf = tf_id,
                 study = 'STES',
                 targets_only = TRUE)

cor_tf <- mutate(cor_tf, cor = cor * 100) %>%
    dcast(tf + feature ~ study, value.var = 'cor')
targets_tf <- select(cor_tf, tf, feature) %>%
    unique() %>%
    mutate(study = 'STES')

dbDisconnect(db)

db <- dbConnect(SQLite(), 'test2.db')
dbWriteTable(db, 'cor_mir', cor_mir, overwrite = TRUE)
dbWriteTable(db, 'targets_mir', targets_mir, overwrite = TRUE)
dbWriteTable(db, 'cor_tf', cor_tf, overwrite = TRUE)
dbWriteTable(db, 'targets_tf', targets_tf, overwrite = TRUE)
dbDisconnect(db)

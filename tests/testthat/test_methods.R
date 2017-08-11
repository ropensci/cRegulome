# load libraries
library(RCurl)
library(R.utils)
library(cRegulome)
context('Testing methods output')

# get db file and make object
get_db(test = TRUE)
gunzip('cRegulome.db.gz')
dat <- get_mir('hsa-let-7b', 'ACC', min_cor = .3, max_num = 5)
cmir <- cmicroRNA(dat)
tidy_cmir <- tidy(cmir)

test_that('tidy output a data.frame identical to output of get_mir', {
    expect_equal(dat, tidy_cmir)
    })

# clean up
if(file.exists('cRegulome.db')) unlink('cRegulome.db')

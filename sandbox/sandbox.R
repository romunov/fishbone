roxygen2::roxygenize()
devtools::load_all()

library(data.table)
# library(gridExtra)
# devtools::use_package("data.table")

mt <- fread("./data/pars.csv", dec = ",",
            colClasses = list(character = c(1, 2), numeric = c(3, 4, 5, 6)),
            stringsAsFactors = FALSE, header = TRUE)
x <- fread("./sandbox/EM.1CPA.txt",
           colClasses = list(character = c(1, 4, 5, 7)))

load("../ngs_pipelines/DAB/data/genotypes_dab_hiseq2_cleaned.RData")

x <- gt[Sample_Name %in% c("EX.1JF8") & Marker %in% c("03")
   , ][, callAllele(c(.BY, .SD), tbase = mt), by = .(Sample_Name, Marker, Plate)]
x <- x[, 4:ncol(x)]
x
system.time(out <- x[, callAllele(c(.BY, .SD), tbase = mt), by = .(Sample_Name, Marker, Plate)])

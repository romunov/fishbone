roxygen2::roxygenize()
devtools::load_all()

library(data.table)
library(gridExtra)
# devtools::use_package("data.table")

mt <- fread("./data/pars.csv", dec = ",",
            colClasses = list(character = c(1, 2), numeric = c(3, 4, 5, 6)),
            stringsAsFactors = FALSE, header = TRUE)
x <- fread("./sandbox/EM.1CPA.txt",
           colClasses = list(character = c(1, 4, 5, 7), numeric = c(2, 3, 6)))

head(x)
system.time(out <- x[, callAllele(c(.BY, .SD), tbase = mt), by = .(Sample_Name, Marker, Plate)])

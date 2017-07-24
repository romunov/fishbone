roxygen2::roxygenize()
devtools::load_all()

library(data.table)
library(gridExtra)
# devtools::use_package("stringi")

# extract data for one sample
# library(readxl)
# xy <- read_excel("./data/genotypes_dinalpbear_gatc_notrash.xlsx")
# xy$Marker <- sprintf("%02d", xy$Marker)
# xy$Run_Name <- NULL
# xy$Sample_ID <- NULL
# xy$TagCombo <- NULL
# xy$Position <- NULL
#
# m1 <- as.data.frame(xy[(xy$Sample_Name %in% "M0P0K.MM") & (xy$Plate %in% 1:8), ])
# m1$Run_Name <- NULL
# m1$Sample_ID <- NULL
# m1$TagCombo <- NULL
# m1$old_Allele <- NULL
# m1$Position <- NULL
# x <- split(m1, f = list(m1$Marker))
# names(x) <- paste("m.", names(x), sep = "")
# me <- list2env(x)
# save(me, file = "./sandbox/samples.RData")

mt <- fread("./data/pars.csv", dec = ",",
            colClasses = list(character = c(1, 2), numeric = c(3, 4, 5, 6)),
            stringsAsFactors = FALSE, header = TRUE)
x <- fread("./sandbox/EM.1CPA.txt",
           colClasses = list(character = c(1, 4, 5, 7), numeric = c(2, 3, 6)))

x <- split(x, f = list(x$Sample_Name, x$Marker, x$Plate))

callAllele(fb = x[[1]], tbase = mt)
callAllele(fb = x[[14]], tbase = mt)
callAllele(fb = x[[27]], tbase = mt)
callAllele(fb = x[[40]], tbase = mt)

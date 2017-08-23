roxygen2::roxygenize()
devtools::load_all()

library(data.table)
# library(gridExtra)
# devtools::use_package("data.table")

mt <- fread("./data/pars.txt", dec = ",",
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

"M0TCU"
outsmp1 <- callAllele(gt[Sample_Name %in% c("M0X0E") & Plate == 3 & Marker == "03", ], tbase = mt)

smp_m0tcu <- gt[Sample_Name %in% c("M0X0E"), ]
outsmp1 <- smp_m0tcu[, callAllele(c(.BY, .SD), tbase = mt), by = .(Sample_Name, Marker, Plate)]
outsmp1 <- outsmp1[, 4:ncol(outsmp1)]
save(smp_m0tcu, file = "./data/smp_m0tcu.RData")
write.table(outsmp1, "test.txt", sep = "\t", row.names = FALSE)

save(mt, file = "./data/mt.RData")

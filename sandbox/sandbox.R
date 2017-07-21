roxygen2::roxygenize()
devtools::load_all()

library(data.table)
library(gridExtra)
# devtools::use_package("gridExtra")

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

tbase <- read.table("sandbox/thresholds.csv", sep = "\t", dec = ",", header = TRUE,
                    colClasses = c("character", "character", "numeric"))
mymotif <- read.table("sandbox/motifs.csv", sep = "\t", header = TRUE,
                      colClasses = "character")
#
load("./sandbox/samples.RData")
x <- as.fishbone(me$m.03.8, motif = mymotif)
plot(x)

#####
mt <- fread("./data/pars.csv", dec = ",",
            colClasses = list(character = c(1, 2, 3), numeric = c(4, 5, 6, 7)),
            stringsAsFactors = FALSE)
x <- fread("./sandbox/EM.1CPA.txt",
           colClasses = list(character = c(1, 4, 5, 6, 7, 8, 9, 10), numeric = c(2, 3)))

x <- split(x, f = list(x$Sample_Name, x$Marker))

i <- 1
g1 <- x[[i]]
g1 <- split(g1, f = g1$Plate)
g1 <- sapply(g1, as.fishbone, motif = mt, simplify = FALSE)
plot(g1[[1]])
plot(prepareRuns(g1))


g1[[1]]
g2 <- g1[[1]]
g2$Read_Count/max(g2$Read_Count)

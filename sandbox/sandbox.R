roxygen2::roxygenize()
devtools::load_all()
# devtools::use_package("gridExtra")

# extract data for one sample
# library(readxl)
# xy <- read_excel("./data/genotypes_dinalpbear_gatc_notrash.xlsx")
# xy$Marker <- sprintf("%02d", xy$Marker)
#
# m1 <- as.data.frame(xy[(xy$Sample_Name %in% "M0TJT.Q") & (xy$Plate %in% 1:8), ])
# m1$Run_Name <- NULL
# m1$Sample_ID <- NULL
# m1$TagCombo <- NULL
# m1$old_Allele <- NULL
# m1$Position <- NULL
# x <- split(m1, f = list(m1$Marker))
# names(x) <- paste("m.", names(x), sep = "")
# me <- list2env(x)
# save(me, file = "./sandbox/samples.RData")

load("./sandbox/samples.RData")
x <- as.fishbone(me$m.03.1, motif = mymotif)
attributes(x)
plot(x)

tbase <- read.table("sandbox/thresholds.csv", sep = "\t", dec = ",", header = TRUE,
                    colClasses = c("character", "character", "numeric"))
mymotif <- read.table("sandbox/motifs.csv", sep = "\t", header = TRUE,
                      colClasses = "character")

m <- mget(names(me)[grepl("m\\.03", names(me))], envir = me)
m <- sapply(m, FUN = as.fishbone, motif = mymotif, simplify = FALSE)
plotRuns(m)

lapply(m, FUN = function(i) {
  ix <- as.fishbone(i, motif.length = nchar(mymotif[mymotif$locus == unique(i$Marker), "motif"]))
  callGenotype(ix, tbase = tbase, motif = mymotif)
})




callGenotype(fb = x, tbase = tbase, motif = mymotif)

markers <- c("m.03", "m.06", "m.14", "m.16", "m.17", "m.25", "m.51",
             "m.57", "m.63", "m.64", "m.65", "m.67", "m.68")

pdf("./sandbox/all_markers.pdf", onefile = TRUE, width = 10, height = 8)
for (i in markers) {
  message(sprintf("Plotting %s", i))
  xy <- mget(ls(pattern = i, envir = me), envir = me)
  xy <- sapply(xy, as.fishbone, simplify = FALSE, motif.length = ifelse(i == "m.14", 5, 4))
  plotRuns(xy)
}
dev.off()


# poglej kako se obnašajo posamezni lokusi
library(readxl)
xy <- as.data.frame(read_excel("./data/genotypes_dinalpbear_gatc_notrash.xlsx"))
xy$Marker <- sprintf("%02d", xy$Marker)
xy$old_Allele <- NULL
xy$Sample_ID <- NULL
xy$Run_Name <- NULL
xy$TagCombo <- NULL
xy$Position <- NULL

tbase <- read.table("sandbox/thresholds.csv", sep = "\t", dec = ",", header = TRUE,
                    colClasses = c("character", "character", "numeric"))
mymotif <- read.table("sandbox/motifs.csv", sep = "\t", header = TRUE,
                      colClasses = "character")

pdf("all_samples_03.pdf", width = 8, height = 6)
x <- by(data = xy, INDICES = list(xy$Marker, xy$Sample), FUN = function(i) {
  message(sprintf("Plotting %s for sample %s", unique(i$Marker), unique(i$Sample)))
  as.fb <- sapply(split(i, f = i$Plate), as.fishbone, motif = mymotif, simplify = FALSE)
  plotRuns(as.fb)
})
dev.off()

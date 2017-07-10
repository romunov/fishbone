roxygen2::roxygenize()
devtools::load_all()
# devtools::use_package("gridExtra")
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
x <- as.fishbone(me$m.14.1, motif.length = 5)
plot(x)

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

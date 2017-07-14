roxygen2::roxygenize()
devtools::load_all()
# devtools::use_package("gridExtra")

# extract data for one sample
library(readxl)
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

tbase <- data.frame(stat = c("B", "IT", "LRT"))
tbase$L <- "14"
tbase$value <- c(1/3, 0.1, 100)

mymotif <- data.frame(locus = c("03", "06", "14", "16", "17", "25", "51", "57", "63", "64",
                              "65", "67", "68"),
                    motif = c("ctat", "aagg", "tttta", "cttt", "cttt", "cttt", "cttt",
                              "catt", "tcca", "ttta", "gata", "attt", "atct"),
                    stringsAsFactors = FALSE)

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



mot <- "tttta"
st <- x$Sequence[6]
al <- x$Sequence[4]
al2 <- x$Sequence[5]
fst <- sub(mot, "", al)
fst2 <- sub(mot, "", al2)

rbind(strsplit(fst, "")[[1]], strsplit(st, "")[[1]])
sum(strsplit(fst, "")[[1]] != strsplit(st, "")[[1]])

rbind(strsplit(fst2, "")[[1]], strsplit(st, "")[[1]])
sum(strsplit(fst2, "")[[1]] != strsplit(st, "")[[1]])

lengths(gregexpr(mot, al))

out <- xy[xy$Marker == "14", ]
e <- regmatches(x = out$Sequence, m = gregexpr(mot, out$Sequence))
elen <- lengths(e)
table(elen)
as.data.frame(out[which(elen == 5), ])

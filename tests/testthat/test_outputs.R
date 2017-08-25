context("Check inputs")
data(smp1)
data(mt)

test_that("Attached data appears OK", {
  expect_equal(nrow(smp1), 536)
  expect_equal(ncol(smp1), 9)
  expect_equal(unique(smp1$Sample_Name), "M0X0E")
  expect_length(unique(smp1$Plate), 4)
  expect_is(smp1, "data.table")
})

context("Check allele calling")

fb03 <- smp1[smp1$Plate == "3" & smp1$Marker == "03", ]
fb06 <- smp1[smp1$Plate == "3" & smp1$Marker == "06", ]
out03 <- callAllele(fb03, tbase = mt)
out06 <- callAllele(fb06, tbase = mt)

# check correct flag

outall <- smp1[, callAllele(c(.BY, .SD), tbase = mt), by = .(Sample_Name, Marker, Plate)]
outallclean <- outall[, 4:ncol(outall)]

test_that("Is callAllele performin OK", {
  expect_is(out03, "data.table")
  expect_equal(nrow(out03), 3)
  expect_named(out03, c("Sample_Name", "Plate", "Read_Count", "Marker", "Run_Name",
                               "length", "Position", "called", "flag", "stutter", "Sequence",
                               "TagCombo"))
  expect_named(outall, c("Sample_Name", "Marker", "Plate", "Sample_Name", "Plate", "Read_Count",
                         "Marker", "Run_Name", "length", "Position", "called", "flag",
                         "stutter", "Sequence", "TagCombo"))
  expect_named(outallclean, c("Sample_Name", "Plate", "Read_Count", "Marker", "Run_Name",
                              "length", "Position", "called", "flag", "stutter", "Sequence",
                              "TagCombo"))
  expect_equal(nrow(outallclean[called == TRUE]), 90)
  expect_equal(nrow(outallclean[stutter == TRUE]), 90)
  expect_equal(out06[out06$Read_Count == 96, flag], "L")
})

fb03[, fbid := 1:.N, by = .(Sample_Name, Plate, Marker)]
st03 <- fishbone:::findStutter(x = fb03[fb03$length == 63 & fb03$Read_Count > 100, ],
            fb = fb03,
            motif = as.character(mt[Marker == "03", "Repeat"]))

test_that("Check if detecting stutter works", {
  expect_length(st03, 1)
  expect_equal(st03, 1)
})

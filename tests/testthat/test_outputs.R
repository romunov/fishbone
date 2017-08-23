library(fishbone)

context("Check inputs")
data(smp_m0tcu)
data(mt)


test_that("Attached data appears OK", {
  expect_equal(nrow(smp_m0tcu), 536)
  expect_equal(ncol(smp_m0tcu), 9)
  expect_equal(unique(smp_m0tcu$Sample_Name), "M0X0E")
  expect_length(unique(smp_m0tcu$Plate), 4)
  expect_is(smp_m0tcu, "data.table")
})

out03 <- callAllele(smp_m0tcu[smp_m0tcu$Plate == "3" & smp_m0tcu$Marker == "03", ], tbase = mt)
test_that("Is callAllele performin OK", {
  expect_is(out03, "data.table")
  expect_equal(nrow(out03), 3)
  expect_named(out03, c("Sample_Name", "Plate", "Read_Count", "Marker", "Run_Name",
                               "length", "Position", "called", "flag", "stutter", "Sequence",
                               "TagCombo"))
})
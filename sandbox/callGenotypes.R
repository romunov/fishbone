#!/usr/bin/env Rscript

# This script will search for specific files on hard drive, read them and call genotypes.

#### Install packages if missing. ####
needs.pkg <- c("optparse", "data.table")
inst.pkg <- needs.pkg[!needs.pkg %in% installed.packages()]

if (length(inst.pkg) > 0) {
  install.packages(inst.pkg)
}
#####

library(optparse)
option_list <- list(make_option(opt_str = c("-d", "--data"), type = "character", default = NULL,
                                help = "Relative or absolute path to raw data from which genotypes
                                are to be called."),
                    make_option(opt_str = c("-p", "--pars"), type = "character", default = NULL,
                                help = "Relative or absolute path to data holding values on thresholds
                                for all loci being processed."),
                    make_option(opt_str = c("-o", "--output"), type = "character", default = NULL,
                                help = "Relative or absolute path to a file where data is to be stored.
                                Do not forget to specify filename extensions."))

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$data)) stop("Please specify --data argument.")
if (is.null(opt$pars)) stop("Please specify --pars argument.")
if (is.null(opt$output)) stop("Please specify --output argument.")

if (!file.exists(opt$data)) stop(sprintf("%s could not be found.", opt$data))
if (!file.exists(opt$pars)) stop(sprintf("%s could not be found.", opt$pars))

library(data.table)

# Read in data.
xy <- fread(opt$data)
pars <- fread(opt$pars)

## prepare motif table
motif <- pars[, c("Alternative_ID", "Repeat")]
names(motif) <- c("locus", "motif")
motif$motif <- tolower(motif$motif)

# Split it so that we process it by run.
xy <- split(xy, f = c(xy$Sample_Name, xy$Marker, xy$Plate))

xy <- sapply(xy, FUN = as.fishbone, motif = pars)

# Convert to fishbone objects.
# Call genotypes on these objects.
# Format output for writing.
# Dump to file.
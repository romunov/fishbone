---
title: 'fishbone: calling alleles from NGS data'
tags:
  - ngs
  - genetics
  - allele
  - sequence
  - DNA
  - tags
  - primers
  - conservation genetics
authors:
 - name: Roman Luštrik
   orcid: 0000-0002-1410-8253
   affiliation: 1
 - name: Tomaž Skrbinšek
   orcid: 0000-0003-4435-7477
   affiliation: 1
affiliations:
 - name: Biotechnical faculty, University of Ljubljana, Jamnikarjeva 101, SI-1000 Ljubljana, Slovenia
   index: 1
date: October 2018
bibliography: bibliography.bib
---

# Summary
Conservation genetics have progressed by leaps and bounds with the introduction of non invasive sampling (sometimes abbreviated as NGS), i.e. by extracting minute quantities of DNA from cells left behind by animals in question [@waits_noninvasive_2005]. In the past, classical sequencing of microsatellite DNA required scientists to read copious amounts of density graphs and in painstaking manual labor call alleles and construct genotypes.

Methods improved in several ways with the advent of next generation sequencing (also known as NGS), also referred to as high throughput sequencing (HTS). While fragment length is known, more importantly so is the actual microsatellite sequence. This enables scientists to construct genotypes on an even finer scale. Nature of the data is that it comes in text form, which is more scalable than looking at the data purely manually.

Our study involved processing around 4500 samples from brown bears _Ursus arctos_ (unpublished). After running the data through the pipeline prepared by @de_barba_high-throughput_2016, we ended up with up to a million lines of data and little time to process. This package helped us solve the calling of alleles, based on which a genotype for each animal was constructed.

# Brief technical overview
When multiplying short repetitive regions of the genome, also called microsatellites, DNA polymerase tends to slip and produce artefacts which have less repeats than the true allele. These "stutters" can be recognized by being a multiple (1x, 2x, 3x...) of a repeat shorter than the allele. This is the pattern used to "call" alleles. It can get more complicated, where stutter of one allele is actually of the same length as another allele. To help with different scenarios, a parameter object can be passed to the function to fine tune allele calling to your experiment. Package help file (`?callAllele`) and repository README should contain all information you need on the algorithm and how to get started.

# References

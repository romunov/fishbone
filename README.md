Package `fishbone` will take a raw data from an NGS (Next Generation Sequencing) pipeline (de Barba et al., 2016) and output a called allele/stutter output. Basedon these outputs, a genotype can be constructed. This problem will be tackled elsewhere.

Input might look something like this. Below is data for sample `M0X0E`, marker `03` and PCR reaction termed `8`. `Sequence` is the raw sequence of the read(s) and `Read_Count` is the number of times this sequence has been detected for this particular combination of `sample * plate * marker`. `TagCombo` and `Position` are related, where former is the tag combination deposited into the position on a 96 well plate.

```
library(fishbone)

data(mt)  # parameters used for filtering
data(smp1)  # some sample data

smp1[Plate == 6 & Marker == "06", ]

    Sample_Name Plate Read_Count Marker Run_Name length Position
 1:       M0X0E     6         76     06    DAB22     67      001
 2:       M0X0E     6         58     06    DAB22     83      001
 3:       M0X0E     6         26     06    DAB22     75      001
 4:       M0X0E     6         10     06    DAB22     79      001
 5:       M0X0E     6          7     06    DAB22     47      001
 6:       M0X0E     6          6     06    DAB22     63      001
 7:       M0X0E     6          6     06    DAB22     71      001
 8:       M0X0E     6          6     06    DAB22     70      001
 9:       M0X0E     6          6     06    DAB22     71      001
10:       M0X0E     6          5     06    DAB22     71      001
11:       M0X0E     6          4     06    DAB22     91      001
12:       M0X0E     6        378     06    DAB22     71      001
13:       M0X0E     6        209     06    DAB22     87      001
                                                                                       Sequence
 1:                         aaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaaagaagacagattgtaa
 2:         aaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaaagaagacagattgtaa
 3:                 aaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaaagaagacagattgtaa
 4:             aaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaaagaagacagattgtaa
 5:                                             aaggaaggaaggaaggaaggaaggaaggaaaagaagacagattgtaa
 6:                             aaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaaagaagacagattgtaa
 7:                     aaggaaggaaggaaggaaggaaggaaggagggaaggaaggaaggaaggaaggaaaagaagacagattgtaa
 8:                      aggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaaagaagacagattgtaa
 9:                     gaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaaagaagacagattgtaa
10:                     aaggaaggaaggaaggaaggaaggaaggaagggaggaaggaaggaaggaaggaaaagaagacagattgtaa
11: aaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaaagaagacagattgtaa
12:                     aaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaaagaagacagattgtaa
13:     aaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaaagaagacagattgtaa
             TagCombo
 1: caggctaa:acaaccga
 2: caggctaa:acaaccga
 3: caggctaa:acaaccga
 4: caggctaa:acaaccga
 5: caggctaa:acaaccga
 6: caggctaa:acaaccga
 7: caggctaa:acaaccga
 8: caggctaa:acaaccga
 9: caggctaa:acaaccga
10: caggctaa:acaaccga
11: caggctaa:acaaccga
12: caggctaa:acaaccga
13: caggctaa:acaaccga
```

Below is an example of how you can call allele calling function on a specific combination `Sample * Marker * Plate`. Sequences not deemed alleles are removed.

```
data(mt)  # parameters used to filter out (see algorithm description)

callAllele(smp1[Plate == 6 & Marker == "06", ], tbase = mt)
   Sample_Name Plate Read_Count Marker Run_Name length Position called flag stutter
1:       M0X0E     6        209     06    DAB22     87      001   TRUE        FALSE
2:       M0X0E     6         58     06    DAB22     83      001  FALSE         TRUE
3:       M0X0E     6        378     06    DAB22     71      001   TRUE        FALSE
4:       M0X0E     6         76     06    DAB22     67      001   TRUE   DL    TRUE
5:       M0X0E     6          6     06    DAB22     63      001  FALSE         TRUE
                                                                                  Sequence
1: aaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaaagaagacagattgtaa
2:     aaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaaagaagacagattgtaa
3:                 aaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaaagaagacagattgtaa
4:                     aaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaaagaagacagattgtaa
5:                         aaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaaagaagacagattgtaa
            TagCombo
1: caggctaa:acaaccga
2: caggctaa:acaaccga
3: caggctaa:acaaccga
4: caggctaa:acaaccga
5: caggctaa:acaaccga
```

Allele with count 76 got two flags (see below for possible flags) - `L` for low number of counts (as defined in the parameter file) and `D` because it is considered in disequilibrium with the highest allele (378 of length 71). Notice that this is an artefact because allele of length 63 is a stutter of allele of length 67 which is in turn stutter of allele of length 71. Stutter of a stutter then. This can be filtered out with appropriate parameters, which should be tailored according to your data.

To run allele calling on the entire dataset, split it by whichever combination you want and call the function.

```
xy <- split(smp1, f = list(smp1$Sample_Name, smp1$Marker, smp1$Plate, smp1$Run_Name))
out <- lapply(xy, FUN = callAllele, tbase = mt, verbose = FALSE)
out <- do.call(rbind, out)

Sample_Name Plate Read_Count Marker Run_Name length Position called flag stutter
1:       M0X0E     3        233     03    DAB22     63      001   TRUE        FALSE
2:       M0X0E     3        284     03    DAB22     59      001   TRUE         TRUE
3:       M0X0E     3         27     03    DAB22     55      001  FALSE         TRUE
4:       M0X0E     3         96     06    DAB22     87      001   TRUE    L   FALSE
5:       M0X0E     3         29     06    DAB22     83      001  FALSE         TRUE
---                                                                                 
158:       M0X0E     6         31     65    DAB22     89      001  FALSE         TRUE
159:       M0X0E     6       2036     67    DAB22     89      001   TRUE        FALSE
160:       M0X0E     6         96     67    DAB22     85      001  FALSE         TRUE
161:       M0X0E     6        489     68    DAB22     74      001   TRUE        FALSE
162:       M0X0E     6         56     68    DAB22     70      001  FALSE         TRUE
                                                                                  Sequence
1:                           aaatcctgtaacaaatctatctatctatctatctatctatctatctatctatctatctatctc
2:                               aaatcctgtaacaaatctatctatctatctatctatctatctatctatctatctatctc
3:                                   aaatcctgtaacaaatctatctatctatctatctatctatctatctatctatctc
4:   aaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaaagaagacagattgtaa
5:       aaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaggaaaagaagacagattgtaa
---                                                                                          
158: gaagcaacagggtatagatatatagagatagatagatagatagatagatagatagatagataaagagatttattataaggaattggctc
159: taattttaattttttaaaatttattttaaaatatttatttatttatttatttatttatttatttatttatttatttttgctcaccacac
160:     taattttaattttttaaaatttattttaaaatatttatttatttatttatttatttatttatttatttatttttgctcaccacac
161:                aacttaccaacaaactaatctatctatctatctatctatctatctatctatctatctatctatctacatatata
162:                    aacttaccaacaaactaatctatctatctatctatctatctatctatctatctatctatctacatatata
          TagCombo
1: acacacac:ttacgcca
2: acacacac:ttacgcca
3: acacacac:ttacgcca
4: acacacac:ttacgcca
5: acacacac:ttacgcca
---                  
158: caggctaa:acaaccga
159: caggctaa:acaaccga
160: caggctaa:acaaccga
161: caggctaa:acaaccga
162: caggctaa:acaaccga
```

If you want to use `data.table` syntax, note that a set of columns will be duplicated.

```
system.time(out <- smp1[, callAllele(c(.BY, .SD), tbase = mt, verbose = TRUE),
                        by = .(Sample_Name, Marker, Plate, Run_Name)])
```

One duplicate can easily be removed (not shown).

### Algorithm
Algorithm used to extract alleles from junk is as follows:

```
0. find max allele height
1. if allele has number of reads < L, flag it as "L"
2. see if allele has stutter
2a. if yes, mark as called
2aa. if A in disbalance (A < D), flag as "D"
2ab. mark stutter as such $stutter = TRUE
2b. if no, check AlleleWithNoStutterHeight
2ba. if x > AlleleWithNoStutterHeight, add flag "N"
2bb. if x < AlleleWithNoStutterHeight, ignore allele
3. if number of unflagged alleles is more than 2 (those marked with D are not counted), add flag "M" to all
```

### Supported flags

* L = low amplification threshold (if for some reason, number of total reads is very low, alleles get a flag)
* N = no stutter (if there was enough reads but no stutter was present)
* D = disbalance - alleles not in balance (expecting 1:1 for heterozygotes, those out of balance flagged)
* M = multiple alleles (self explanatory)

### Filtering parameters

`mt` is a parameter data.frame where stutter heights, absolute read count etc. are specified. Notice the addition of extra columns `called`, `flag` and `stutter`. This will depend on your experiment so modify or create a new one as needed.

```
> head(mt)
   Marker Repeat Stutter Disbalance LowCount AlleleWithNoStutterHeight
1:     03   ctat    0.18       0.33      100                       100
2:     06   aagg    0.18       0.33      100                       100
3:     14  tttta    0.18       0.33      100                       100
4:     16   cttt    0.18       0.33      100                       100
5:     17   cttt    0.18       0.33      100                       100
6:     25   cttt    0.18       0.33      100                       100
```

### References
De Barba M, Miquel C, Lobréaux S, et al (2016) High-throughput microsatellite genotyping in ecology: improved accuracy, efficiency, standardization and success with low-quantity and degraded DNA. Molecular Ecology Resources 1–16. doi: 10.1111/1755-0998.12594

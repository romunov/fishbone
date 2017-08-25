Package `fishbone` will take a raw genotype from an NGS (Next Generation Sequencing) pipeline and output a called allele/stutter output.
This output can be further imported into an appropriate application to construct consensus genotypes and statistics.

Input might look something like this. Below is data for sample `M0X0E` for marker `03` and PCR reaction termed `8`. Notice which columns are present. `Sequence` is the raw sequence and `Read_Count` is the number of times this sequence was detected for this particular combination of `sample * plate * marker`. `TagCombo` and `Position` are related, where former is the tag combination deposited into the position on a 96 well plate.

```{r}
library(fishbone)
data(mt)
data(smp1)
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

Below is an example of how you can call allele calling function on a combination `Sample * Marker * Plate`. Sequences not deemed alleles are removed.

```{r}
data(mt) # parameter file used to filter out (see algorithm description)
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
`mt` is a parameter file where stutter heights, absolute read count etc. are specified. Notice the addition of extra columns `called`, `flag` and `stutter`.

Allele with count 76 got two flags (see below for possible flags) - `L` for low number of counts (as defined in the parameter file) and `D` because it is considered in disequilibrium with the highest allele (378 of length 71). Notice that this is an artefact because allele of length 63 is a stutter of allele of length 67 which is in turn stutter of allele of length 71. Stutter of a stutter then. This can be filtered out with appropriate parameters, which should be tailored according to your data.

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

From the above we can deduce that if a sequence also has a corresponding stutter of appropriate height(s), it is called as an allele. Sequences functioning as stutters are marked as such.

Known flags are:

```
* L = low amplification threshold (if for some reason, number of total reads is very low, alleles get a flag)
* N = no stutter (if there was enough reads but no stutter was present)
* D = disbalance - alleles not in balance (expecting 1:1 for heterozygotes, those out of balance flagged)
* M = multiple alleles (self explanatory)
```
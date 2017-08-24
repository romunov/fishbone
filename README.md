Package `fishbone` will take a raw genotype from an NGS (Next Generation Sequencing) pipeline and output a called allele/stutter output.
This output can be further imported into an appropriate application to construct consensus genotypes and statistics.

Input might look something like this:

```{r}
25:     EX.1JF8     8        146     03    DAB16     59      038             aaatcctgtaacaaatctatctatctatctatctatctatctatctatctatctatctc
26:     EX.1JF8     8         77     03    DAB16     63      038         aaatcctgtaacaaatctatctatctatctatctatctatctatctatctatctatctatctc
27:     EX.1JF8     8         15     03    DAB16     55      038                 aaatcctgtaacaaatctatctatctatctatctatctatctatctatctatctc
28:     EX.1JF8     8          5     03    DAB16     59      038             aaatcctgtaacaaatctatctatctatctatctatctatctatctatctatctgtctc
```

and when the function is called on this set (Sample * Marker * Plate combination), junk sequences are removed.

```{r}
x[, callAllele(c(.BY, .SD), tbase = mt), by = .(Sample_Name, Marker, Plate)]

    Sample_Name Plate Read_Count Marker Run_Name length Position called flag stutter
12:     EX.1JF8     8         77     03    DAB16     63      038   TRUE        FALSE     aaatcctgtaacaaatctatctatctatctatctatctatctatctatctatctatctatctc
13:     EX.1JF8     8        146     03    DAB16     59      038   TRUE         TRUE         aaatcctgtaacaaatctatctatctatctatctatctatctatctatctatctatctc
14:     EX.1JF8     8         15     03    DAB16     55      038  FALSE         TRUE             aaatcctgtaacaaatctatctatctatctatctatctatctatctatctatctc
```

`mt` is a parameter file where stutter heights, absolute read count etc. are specified. Notice the addition of extra columns `called`, `flag` and `stutter`.

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

From the above we can deduce that if a sequence also has a corresponding stutter of appropriate height(s), it is called as an allele. Sequences functioning as stutters are marked as such. If there is something of interest in the sequence (count), it receives a flag. Flags can be

```
* L = low amplification threshold (if for some reason, number of total reads is very low, alleles get a flag)
* N = no stutter (if there was enough reads but no stutter was present)
* D = disbalance - alleles not in balance (expecting 1:1 for heterozygotes, those out of balance flagged)
* M = multiple alleles (self explanatory)
```
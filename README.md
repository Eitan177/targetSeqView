% targetSeqView vignette 
% Eitan Halper-Stromberg
% 2013/08/27





# Installation
install directly from github

```r
library(devtools)
install_github("targetSeqView", "Eitan177")
```


# Abstract

This package is designed to evaluate structural rearrangment calls from a candidate list, the output for tools such as HYDRA (Quilan, 2010), GASV (Sindi, 2009), VariationHunter (Hormozdiari, 2010), etc. The user should have a text file with one row per candidate structural rearrangment. For each
candidate rearrangement, read-pairs from the two loci will be read in from a bam file and realigned three different ways.
One of these realignments supports the structural variant, with read-pairs realigned to a sequence representing the rearranged
sequence (the sequence of the two loci concatenated together). The other two realignments support no structural rearrangement,
with read-pairs realigned to the two sequences representing contiguous fragments of the reference genome taken from each of the two loci.

# Quick Start Example

We start with a text file containing 20 candidate deletion junctions and the bam file containing our read alignments from a whole
genome sequencing experiment (for the purposes of this vignette the bam file contains only reads aligning within the regions that we will
interrogate). The text file contains the loci allegedely involved in each deletion, one deletion per row. For each row we will load
reads aligning to the two loci involved in the alleged deletion from our bam file, realign these reads, and calculate our
likelihood score.



```r
suppressMessages(library(targetSeqView))
path <- system.file("extdata", package = "targetSeqView")
## This method utilizes the foreach package for parallelization, set nodes to
## however many cpus are available.
nodes = 1
registerDoMC(nodes)
## create an instance of the candidates class
candidateDels <- new("candidates")
## set the path where bam files are located (if not in the currect working
## directory)
bamFilePath(candidateDels) <- path
## set the name of the text file containing candidate SVs (full path if not
## in the working directory)
candidatesFileName(candidateDels) <- file.path(path, "wholeGenomeDeletionCandidates.txt")
## set the build of the (human) genome
build(candidateDels) <- "hg19"
## set the read length
readLength(candidateDels) <- 101
## set the mismatch rate for each position along the read length
mmRate(candidateDels) <- precomputedWholeGenome101bpMMRate()
## set the indel rate for reach position along the read length
indelRate(candidateDels) <- precomputedWholeGenome101bpIndelRate()
```


note: mismatch and indel rates may be calculated based upon reads from a bam file containing
normal alignments, the bamFile argument should contain the full path with the bam name if the
bam is not in the current directory. The following 3 lines are unevaluated

```r
normalBam <- "Path/To/Normal/bamfile.bam"
errorRates <- getErrorRate(normalBam)
mmRate(candidateDels) <- errorRates[["mmRate"]]
indelRate(candidateDels) <- errorRates[["indelRate"]]
```


We first obtain likelihood scores for candidates without performing full smith-waterman
realignment on all reads for all 3 alignement configurations. We instead use alignment
information in the cigar strings and md tags in our bam file to obtain mismatches and indels
for the alignments supporting SVs. In addition, We forgo, for the moment, returning a
data.frame formatted for our plot function. This should take a few (1-5) seconds per candidate

```r
candidateDels <- quickScore(candidateDels, verbose = TRUE)
```

```
## Of 20 events now working on  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
```

```r
## view values returned
print(candidateDels@quickScore)
```

```
##        1        2        3        4        5        6        7        8 
## 1280.592  896.235  881.604  793.243  556.465  844.628  741.042  461.249 
##        9       10       11       12       13       14       15       16 
## 1065.448  749.678 -252.563  -77.755  -41.597  -35.835  -33.225  -32.736 
##       17       18       19       20 
##  -28.493  -25.971  -25.165   -7.992
```

```r
### In this case we have validation data for these candidates
indexOfvalidated <- 1:10
validated <- candidateDels@quickScore[indexOfvalidated]
failedvalidation <- candidateDels@quickScore[-indexOfvalidated]
```


### Figure 1: boxplot of likelihood scores from 20 candidates; 10 validated, 10 failed validation

```r
boxplot(list(validated = validated, failed = failedvalidation, all = candidateDels@quickScore), 
    ylab = "log likelihood score")
```

![Distribution of 20 candidate deletions taken from a whole-genome sequencing dataset, broken down by validation status](figure/graphics.png) 


# Scoring and Viewing

In this section we will view read alignments at the junctions of 3 candidate structural variants. We will use a sequencing dataset
taken from a target-capture exerpiment. As with the last section, we start with a text file containing the candidate structural
variants (in this case 1 inversion and 2 chromosomal translocations) and a bam file containing read alignments.

```r
## create an instance of the candidates class
candidateSVs <- new("candidates")
bamFilePath(candidateSVs) <- path
candidatesFileName(candidateSVs) <- file.path(path, "targetCaptureSVs.txt")
build(candidateSVs) <- "hg19"
readLength(candidateSVs) <- 100
mmRate(candidateSVs) <- precomputedTargetCapture100bpMMRate()
indelRate(candidateSVs) <- precomputedTargetCapture100bpIndelRate()

## fullScoreAndview will perform full smith-waterman realignment for all
## reads in the 3 configurations. In addition, if the input text file
## contains a SplitsSample column, the function will look for split-reads
## within the bam file specified by the column 'SplitsSample'
candidateSVs <- fullScoreAndView(candidateSVs, verbose = TRUE, findSplitReads = TRUE)
```

```
## [1] "Working on event 1 of 3"
## [1] "primary alignment for event 1 done"
## [1] "secondary alignment (1 of 2) for event 1 done"
## [1] "secondary alignment (2 of 2) for event 1 done"
## [1] "Working on event 2 of 3"
## [1] "primary alignment for event 2 done"
## [1] "secondary alignment (1 of 2) for event 2 done"
## [1] "secondary alignment (2 of 2) for event 2 done"
## [1] "Working on event 3 of 3"
## [1] "primary alignment for event 3 done"
## [1] "secondary alignment (1 of 2) for event 3 done"
## [1] "secondary alignment (2 of 2) for event 3 done"
```

```r
## The Scores for these 3 candidates:
print(candidateSVs@fullScore)
```

```
## [1]   1.182 957.485 501.009
```


Let's view the negative. Read-pair alignments supporting the SV look good, but read-pair
alignments supporting contiguous fragments also look good, especially for side 1. This is how we identify the event as a negative and explains why the event recieved a low likelihood score 

### Figure 2: A chromosomal translocation that failed to validate

```r
plotSV(candidateSVs, indices = 1, flipLeftandRight = TRUE)
```

![A negative (i.e not real) chromosomal translocation. The top plot shows read-pair alignments supporting the SV and the bottom plots show read-pair alignments supporting contiguous sequences.](figure/graphics2.png) 

```
## .
```


Let's view the first positive. Read-pair alignments supporting the SV look good, read-pair
alignments supporting contiguous fragments do not look good because in both contiguous fragment
alignment pictures, one read from each pair has many mismatches/indels

### Figure 3: A validated inversion

```r
plotSV(candidateSVs, indices = 2)
```

![A positive (i.e real) inversion. The top plot shows read-pair alignments supporting the SV and the bottom plots show read-pair alignments supporting contiguous sequences](figure/graphics3.png) 

```
## .
```


Lets view the second positive. Read-pair alignments supporting the SV look good, albeit for a
different reason than the first positive. In this picture we have some reads aligning well to the
chr14 side and their partners aligning across the junction of the SV (i.e split-reads). The
split-reads align well to both sides. There are a few mismatches right at the junction for these
split reads but otherwise they match the reference. The contiguous fragment alignments do not look
good, as we would expect. Again, the flipLeftandRight option is a style preference, putting the
chr14 junction on the left and the chr18 junction on the right.

### Figure 4: A validation chromosomal translocation

```r
plotSV(candidateSVs, indices = 3, flipLeftandRight = TRUE)
```

![A positive (i.e real) chromosomal translocation. The top plot shows read-pair alignments supporting the SV and the bottom plots show read-pair alignments supporting contiguous sequences.](figure/graphics4.png) 

```
## .
```



```r
print(sessionInfo(), locale = FALSE)
```

```
## R version 3.0.1 (2013-05-16)
## Platform: x86_64-apple-darwin10.8.0 (64-bit)
## 
## attached base packages:
## [1] grid      parallel  stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] BSgenome.Hsapiens.UCSC.hg19_1.3.19 targetSeqView_0.99                
##  [3] doMC_1.3.0                         iterators_1.0.6                   
##  [5] ggplot2_0.9.3.1.99                 BSgenome_1.28.0                   
##  [7] Rsamtools_1.12.4                   Biostrings_2.28.0                 
##  [9] GenomicRanges_1.12.5               IRanges_1.18.3                    
## [11] BiocGenerics_0.6.0                 foreach_1.4.1                     
## [13] knitr_1.4.1                       
## 
## loaded via a namespace (and not attached):
##  [1] bitops_1.0-6       codetools_0.2-8    colorspace_1.2-2  
##  [4] compiler_3.0.1     dichromat_2.0-0    digest_0.6.3      
##  [7] evaluate_0.4.7     formatR_0.9        gtable_0.1.2      
## [10] labeling_0.2       MASS_7.3-28        munsell_0.4.2     
## [13] plyr_1.8           proto_0.3-10       RColorBrewer_1.0-5
## [16] reshape2_1.2.2     scales_0.2.3       stats4_3.0.1      
## [19] stringr_0.6.2      tools_3.0.1        zlibbioc_1.6.0
```

# Bibliography
- Hormozdiari F, Hajirasouliha I, Dao P, Hach F, Yorukoglu D, Alkan C, Eichler E
and Sahinalp S (2010). “Next-Generation Variationhunter: Combinatorial Algorithms
For Transposon Insertion Discovery.” _Bioinformatics_, *26*, pp. i350-i357. ISSN
1367-4803, <URL: http://dx.doi.org/10.1093/bioinformatics/btq216>.
- Quinlan A, Clark R, Sokolova S, Leibowitz M, Zhang Y, Hurles M, Mell J and Hall I
(2010). “Genome-Wide Mapping And Assembly of Structural Variant Breakpoints in
The Mouse Genome.” _Genome Research_, *20*, pp. 623-635. ISSN 1088-9051, <URL:
http://dx.doi.org/10.1101/gr.102970.109>.
- Sindi S, Helman E, Bashir A and Raphael B (2009). “A Geometric Approach For
Classification And Comparison of Structural Variants.” _Bioinformatics_, *25*,
pp. i222-i230. ISSN 1367-4803, <URL:
http://dx.doi.org/10.1093/bioinformatics/btp208>.

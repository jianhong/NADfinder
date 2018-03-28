# CHANGES IN VERSION 1.3.5

* Speed up tileCount per Herve's recommendation
* Add Herve Pages as a coauthor
* resave triplicate.count and single.count
 
# CHANGES IN VERSION 1.3.3

* Fix the warning in windows for the help files.
* Add Haibo Liu as a coauthor.
* Add different methods for transformation and combining p-values.

# CHANGES IN VERSION 1.3.1

* remove Bessel's correction for zscore calculation
* narrow the distance when trim the peaks to avoid overlaps.
* change the baseline.modpolyfit from degree 4 to degree 3.

# CHANGES IN VERSION 1.1.2

* update the documentations.

# CHANGES IN VERSION 1.1.1

* add parameter in backgroundCorrection.

# CHANGES IN VERSION 0.99.4

* replace slidingWindows by GenomicRanges::slidingWindows

# CHANGES IN VERSION 0.99.3

* add new function cumulativePercentage
* add DNASeq, PeakDetection to biocViews.
* update documentations of vignette by provide a narrative better explaining 
use of the package.
* indicate argument type for each parameter in the documentations.
* use vectorized rather than iteractive operations in place of other uses of
`apply()`
* use `seq_len` and `seq_along` to replace `:`

# CHANGES IN VERSION 0.99.2

* change GRanges to RangedSummarizedExperiment
* update documentations
* change plotCorrelations to getCorrelations
* change sig2bedgraph to exportSignals
* change tileGRanges to slidingWindows

# CHANGES IN VERSION 0.99.1

* change NEWS to NEWS.md
* add details to README.md
* update documentations

# CHANGES IN VERSION 0.99.0

* update the version to 0.99.0

# CHANGES IN VERSION 0.1.4

* add new function plotSig

# CHANGES IN VERSION 0.1.3

* add new function plotCorrelations

# CHANGES IN VERSION 0.1.2

* add new authors.

# CHANGES IN VERSION 0.1.1

* add the power to analyze duplicates.

# CHANGES IN VERSION 0.1.0

* Create the package.

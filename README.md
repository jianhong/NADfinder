# NADfinder

[![platforms](http://bioconductor.org/shields/availability/3.5/NADfinder.svg)](http://bioconductor.org/packages/devel/bioc/html/NADfinder.html)
[![build](http://bioconductor.org/shields/build/devel/bioc/NADfinder.svg)](http://bioconductor.org/packages/devel/bioc/html/NADfinder.html)
[![test coverage](https://codecov.io/github/Bioconductor-mirror/NADfinder/branch/master)](http://bioconductor.org/packages/devel/bioc/html/NADfinder.html)

Call wide peaks for sequencing data

Nucleoli serve as major organizing hubs for the three-dimensional structure 
of mammalian heterochromatin. Specific loci, termed
nucleolar-associated domains (NADs), form frequent three-dimensional 
associations with nucleoli. Early mammalian 
development is a critical period to study NAD biological function, because 
interactions between pericentromeric chromatin and perinucleolar regions are 
particularly dynamic during embryonic development . 
We therefore propose for the first time to map the 
NADs in the mouse genome, determine how these associations are altered during 
embryonic stem cell (ESC) differentiation, and develop tools for study of 
these higher-order chromosome interactions in fixed and live single cells. 

## Installation

To install this package, start R and enter:

```r
library(BiocInstaller)
biocLite("NADfinder")
```

## Documentation

To view documentation of NADfinder, start R and enter:
```r
browseVignettes("NADfinder")
```

## Why NADfinder

The peaks are wide peaks (over 10K bps) and the chromosomes were acrocentric. 

## Steps

1. Reads count: We move the window (w) along the genome with step (s) and count
the reads in each window. This step can smooth the coverage, 
which is a good for wide peaks.
2. Ratio calculation: ratio = log2(nucleosome counts / genome counts), 
pseudocount will be used to avoid x/0 by x/pseudocount.
3. Background correction: Because the the ratios are higher in 5’end than 
3’ end, we applied modified polynomial fitting to remove the background. 
With this step, the baseline of the ratios will keep at 0 along each 
chromosome. More details could refer: CHAD A. LIEBER and ANITA 
MAHADEVAN-JANSEN: Automated Method for Subtraction of 
Fluorescence from Biological Raman Spectra 
(http://journals.sagepub.com/doi/pdf/10.1366/000370203322554518).  
4. Curve smooth: smoothed curve will be used for peak range detection. 
We applied butterworth filter to smooth the ratio curve. 
The idea is that with this filter, the high frequency noise will be removed.
5. Visualization of the background correction and signal processing filter
to double check the process is reasonable.
6. Calculate z-score for each chromosome by smoothed ratios and call peaks. 
Because the peaks are start from previous valley to next valley, peaks will be 
trimmed by background corrected signals from both shoulder of the curve to 
make sure the peak region  does not indclude the parts of valley.   
7. Export the peaks into bigwig files for visualization.
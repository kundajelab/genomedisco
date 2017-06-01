# genomedisco
[![Build Status](https://travis-ci.org/kundajelab/genomedisco.svg?branch=master)](https://travis-ci.org/kundajelab/genomedisco)
[![license](https://img.shields.io/github/license/mashape/apistatus.svg?maxAge=2592000)](https://github.com/kundajelab/genomedisco/blob/master/LICENSE)

`genomedisco` is a package for comparing contact maps of 3D genome structures, obtained from experiments such as Hi-C, Capture-C, ChIA-PET, HiChip, etc. It uses graph diffusion to smooth contact maps, and then compares them, resulting in a reproducibility score that can be used for quality control of biological replicates.

Installation
===

1. Install [Anaconda](https://www.continuum.io/downloads). 
2. Get genomedisco with the following command:
```
git clone http://github.com/kundajelab/genomedisco
```
genomedisco is compatible with Python 2.

Quick start
====

Say you want to compare 2 contact maps. For this example, we will use a subset of datasets from Rao et al., 2014. 

```
genomedisco all --metadata_samples examples/metadata.samples --metadata_pairs examples/metadata.pairs --nodes examples/Nodes.w40000.bed.gz --outdir examples/output 
```

The analysis produces a beautiful html report of the results. For the example we just ran, the html is here: http://htmlpreview.github.io/?http://github.com/kundajelab/genomedisco/blob/master/examples/output/results/sample1.vs.sample2/report.sample1.vs.sample2.genomedisco.html

Running reproducibility analysis in batches
====

In the above example, we computed reproducibility for comparing 2 samples. But what if you have multiple comparisons of interest? 

All you need to do is modify the metadata files (in our example `examples/metadata.batch.samples` and `examples/metadata.batch.pairs`), and then run the similar command:

```
genomedisco all --metadata_samples examples/metadata.batch.samples --metadata_pairs examples/metadata.batch.pairs --nodes examples/Nodes.w40000.bed.gz --outdir examples/output 
```

Again, you can visualize the pretty reports as html files. Note that sample1 and sample2 are from the same cell type, while sample3 is a different cell type. This is reflected in the reproducibility scores. 
- sample1 vs sample2: http://htmlpreview.github.io/?http://github.com/kundajelab/genomedisco/blob/master/examples/output_batch/results/sample1.vs.sample2/report.sample1.vs.sample2.genomedisco.html
- sample1 vs sample3: http://htmlpreview.github.io/?http://github.com/kundajelab/genomedisco/blob/master/examples/output_batch/results/sample1.vs.sample3/report.sample1.vs.sample3.genomedisco.html
- sample2 vs sample3: http://htmlpreview.github.io/?http://github.com/kundajelab/genomedisco/blob/master/examples/output_batch/results/sample2.vs.sample3/report.sample2.vs.sample3.genomedisco.html

More questions?
====
Contact Oana Ursu

oursu@stanford.edu

We are excited to hear from you and how we can improve GenomeDISCO!


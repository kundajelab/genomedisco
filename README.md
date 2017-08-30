# GenomeDISCO

`GenomeDISCO` (DIfferences between Smoothed COntact maps) is a package for comparing contact maps of 3D genome structures, obtained from experiments such as Hi-C, Capture-C, ChIA-PET, HiChip, etc. It uses random walks on the contact map graph for smoothing before comparing the contact maps, resulting in a concordance score that can be used for quality control of biological replicates.

Read the full paper here: 
Oana Ursu, Nathan Boley, Maryna Taranova, Y. X. Rachel Wang, Galip Gurkan Yardimci, William Stafford Noble, Anshul Kundaje. GenomeDISCO: A concordance score for chromosome conformation capture experiments using random walks on contact map graphs. bioRxiv: http://www.biorxiv.org/content/early/2017/08/29/181842

Installation
===

1. Install [Anaconda](https://www.continuum.io/downloads). 
2. Obtain and install GenomeDISCO with the following commands:
```
git clone http://github.com/kundajelab/genomedisco
genomedisco/install_scripts/install_genomedisco.sh
```

GenomeDISCO is compatible with Python 2.

Quick start
====

Say you want to compare 2 contact maps. For this example, we will use a subset of datasets from Rao et al., 2014. 

```
cd genomedisco
#configure example
examples/configure_example.sh

#run analysis
python reproducibility_analysis/chromatin3d_replicateQC.py run_all --metadata_samples examples/metadata.samples --metadata_pairs examples/metadata.pairs --nodes examples/Nodes.w40000.bed.gz --outdir examples/output 
```

The analysis produces a beautiful html report of the results. For the example we just ran, the html is here: http://htmlpreview.github.io/?http://github.com/kundajelab/genomedisco/blob/master/examples/output/results/sample1.vs.sample2/report.sample1.vs.sample2.genomedisco.html

To run reproducibility analysis in batches (more than one comparison), see the documentation.

Documentation
=============

See the full documentation here.

More questions?
====
Contact Oana Ursu

oursu@stanford.edu

We are excited to hear from you and how we can improve GenomeDISCO!


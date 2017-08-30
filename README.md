# GenomeDISCO



`GenomeDISCO` (**DI**fferences between **S**moothed **CO**ntact maps) is a package for comparing contact maps of 3D genome structures, obtained from experiments such as Hi-C, Capture-C, ChIA-PET, HiChip, etc. It uses random walks on the contact map graph for smoothing before comparing the contact maps, resulting in a concordance score that can be used for quality control of biological replicates.

Read the full paper here: 
*A concordance score for chromosome conformation capture experiments using random walks on contact map graphs.* Oana Ursu, Nathan Boley, Maryna Taranova, Y. X. Rachel Wang, Galip Gurkan Yardimci, William Stafford Noble, Anshul Kundaje. bioRxiv: http://www.biorxiv.org/content/early/2017/08/29/181842

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

First, configure the files used in the example:

```
genomedisco/examples/configure_example.sh
```

Then run the concordance analysis:

```
python genomedisco/reproducibility_analysis/chromatin3d_replicateQC.py run_all --method GenomeDISCO --metadata_samples examples/metadata.samples --metadata_pairs examples/metadata.pairs --nodes examples/Nodes.w40000.bed.gz --outdir examples/output 
```

To run reproducibility analysis in batches (more than one comparison), see the documentation.

Documentation
=============

See the full documentation here.

More questions?
====
Contact Oana Ursu

oursu@stanford.edu

We are excited to hear from you and how we can improve GenomeDISCO!


# GenomeDISCO



`GenomeDISCO` (**DI**fferences between **S**moothed **CO**ntact maps) is a package for comparing contact maps of 3D genome structures, obtained from experiments such as Hi-C, Capture-C, ChIA-PET, HiChip, etc. It uses random walks on the contact map graph for smoothing before comparing the contact maps, resulting in a concordance score that can be used for quality control of biological replicates.

Read the full paper here: 
*GenomeDISCO: A concordance score for chromosome conformation capture experiments using random walks on contact map graphs.* Oana Ursu, Nathan Boley, Maryna Taranova, Y. X. Rachel Wang, Galip Gurkan Yardimci, William Stafford Noble, Anshul Kundaje. bioRxiv: http://www.biorxiv.org/content/early/2017/08/29/181842

Installation
===

1. Install [Anaconda](https://www.continuum.io/downloads). GenomeDISCO is compatible with Python 2.
2. Obtain and install GenomeDISCO with the following commands:
```
git clone http://github.com/kundajelab/genomedisco
genomedisco/install_scripts/install_genomedisco.sh
```

**Note if you are installing these locally**: There are a few parameters you can provide to the installation script, to point it to your desired python installation, R installation, R library, modules and bedtools installation. Thus, you can run the above script as follows:

```
genomedisco/install_scripts/install_genomedisco.sh --pathtopython /path/to/your/python --pathtor /path/to/your/R --rlib /path/to/your/Rlibrary --modules modulename --pathtobedtools path/to/your/bedtools
```

Quick start
====

Say you want to compare 2 contact maps. For this example, we will use a subset of datasets from Rao et al., 2014. 

First, configure the files used in the example:

```
genomedisco/examples/configure_example.sh
```

Then run the concordance analysis:

```
cd genomedisco
python reproducibility_analysis/3DChromatin_ReplicateQC.py run_all --method GenomeDISCO --metadata_samples examples/metadata.samples --metadata_pairs examples/metadata.pairs --bins examples/Nodes.w40000.bed.gz --outdir examples/output 
```

For detailed explanations of all inputs to GenomeDISCO, go to the ["Inputs" section of the documentation](#iputs)

To run reproducibility analysis in batches (more than one comparison), all you need to do is modify the `--metadata_samples` and `--metadata_pairs` to add the additional samples and sample pairs respectively that you wish to compare.

Running other methods for measuring concordance and QC of Hi-C data
====

**coming soon**

GenomeDISCO supports computing concordance scores for Hi-C data using not only the GenomeDISCO framework, but also HiCRep (http://github.com/qunhualilab/hicrep), HiC-Spector (http://github.com/gersteinlab/HiC-spector) and QuASAR-Rep (part of the hifive suite at http://github.com/bxlab/hifive). In addition, it also computes QC scores for Hi-C data using QuASAR-QC (part of the hifive suite at http://github.com/bxlab/hifive). Thanks to Tao Yang and Michael Sauria for providing wrapper scripts around their methods.

Documentation
=============

Inputs
-----

See the full documentation here.

More questions?
====
Contact Oana Ursu

oursu@stanford.edu

We are excited to hear from you and how we can improve GenomeDISCO!


# GenomeDISCO



`GenomeDISCO` (**DI**fferences between **S**moothed **CO**ntact maps) is a package for comparing contact maps of 3D genome structures, obtained from experiments such as Hi-C, Capture-C, ChIA-PET, HiChip, etc. It uses random walks on the contact map graph for smoothing before comparing the contact maps, resulting in a concordance score that can be used for quality control of biological replicates.

Read the full paper here: 
*GenomeDISCO: A concordance score for chromosome conformation capture experiments using random walks on contact map graphs.* Oana Ursu, Nathan Boley, Maryna Taranova, Y. X. Rachel Wang, Galip Gurkan Yardimci, William Stafford Noble, Anshul Kundaje. Bioinformatics, 2018. https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/bty164/4938489?redirectedFrom=fulltext


Installation
===

1. Install [Anaconda](https://www.continuum.io/downloads). GenomeDISCO is compatible with Python 2.
2. Obtain and install GenomeDISCO with the following commands:
```
git clone http://github.com/kundajelab/genomedisco
pip install --editable genomedisco/
```

Quick start
====

Say you want to compare 2 contact maps. For this example, we will use a subset of datasets from Rao et al., 2014. 

First, configure the files used in the example: (this will create all input files necessary for the example on which we will run GenomeDISCO)

```
genomedisco/examples/configure_example.sh
```

Then run the concordance analysis:

```
cd genomedisco
genomedisco run_all --metadata_samples examples/metadata.samples --metadata_pairs examples/metadata.pairs --bins examples/Bins.w50000.bed.gz --outdir examples/output 
```

For detailed explanations of all inputs to GenomeDISCO, see the ["Inputs" section below](#inputs)

To run reproducibility analysis in batches (more than one comparison), all you need to do is modify the `--metadata_samples` and `--metadata_pairs` to add the additional samples and sample pairs respectively that you wish to compare. For details, see ["Analyzing multiple dataset pairs"](#analyzing-multiple-dataset-pairs)

Running other methods for measuring concordance and QC of Hi-C data
====

To run other available methods for computing the reproducibility of Hi-C data, refer to the repository http://github.com/kundajelab/3DChromatin_ReplicateQC and follow the instructions there.

The reproducibility methods supported in 3DChromatin_ReplicateQC are:
- GenomeDISCO (http://github.com/kundajelab/genomedisco)
- HiCRep (http://github.com/qunhualilab/hicrep) 
- HiC-Spector (http://github.com/gersteinlab/HiC-spector) 
- QuASAR-Rep (part of the hifive suite at http://github.com/bxlab/hifive) 

Note: given that both GenomeDISCO and 3DChromatin_ReplicateQC use the same underlying base code, they share the parameter options below, resulting in shared README sections for these.

Inputs
=============

Before running GenomeDISCO, make sure to have the following files:

- **contact map** For each of your samples, you need a file containing the counts assigned to each pair of bins in your contact map, and should have the format `chr1 bin1 chr2 bin2 value`. Note: GenomeDISCO assumes that this file contains the contacts for all chromosomes, and will split it into individual files for each chromosome.

- **bins** This file contains the full set of genomic regions associated with your contact maps, in the format `chr start end name` where name is the name of the bin as used in the contact map files above. GenomeDISCO supports both fixed-size bins and variable-sized bins (e.g. obtained by partitioning the genome into restriction fragments). 

GenomeDISCO takes the following inputs:

- `--metadata_samples` Information about the samples being compared. Tab-delimited file, with columns "samplename", "samplefile". Note: each samplename should be unique. Each samplefile listed here should follow the format "chr1 bin1 chr2 bin2 value

- `--metadata_pairs` Each row is a pair of sample names to be compared, in the format "samplename1 samplename2". Important: sample names used here need to correspond to the first column of the --metadata_samples file.

- `--bins` A (gzipped) bed file of the all bins used in the analysis. It should have 4 columns: "chr start end name", where the name of the bin corresponds to the bins used in the contact maps.

- `--re_fragments` Add this flag if the bins are not uniform bins in the genome (e.g. if they are restriction-fragment-based).By default, the code assumes the bins are of uniform length.

- `--parameters_file` File with parameters for reproducibility and QC analysis. For details see ["Parameters file"](#parameters-file)

- `--outdir` Name of output directory. DEFAULT: replicateQC

- `--running_mode` The mode in which to run the analysis. This allows you to choose whether the analysis will be run as is, or submitted as a job through sge or slurm. Available options are: "NA" (default, no jobs are submitted). Coming soon: "sge", "slurm"

- `--concise_analysis` Set this flag to obtain a concise analysis, which means replicateQC is measured but plots that might be more time/memory consuming are not created. This is useful for quick testing or running large-scale analyses on hundreds of comparisons.

- `--subset_chromosomes` Comma-delimited list of chromosomes for which you want to run the analysis. By default the analysis runs on all chromosomes for which there are data. This is useful for quick testing

Analyzing multiple dataset pairs
======
To analyze multiple pairs of contact maps, all you need to do is add any additional datasets you want to analyze to the `--metadata_samples` file and any additional pairs of datasets you want to compare to the `--metadata_pairs` files. 

Parameters file
======

The parameters file specifies the parameters to be used with GenomeDISCO (and any of the other methods GenomeDISCO supports). The format of the file is: `method_name parameter_name parameter_value`. The default parameters file used by GenomeDISCO is:

```
GenomeDISCO|subsampling	lowest
GenomeDISCO|tmin	3
GenomeDISCO|tmax	3
GenomeDISCO|norm	sqrtvc
GenomeDISCO|scoresByStep	no
GenomeDISCO|removeDiag	yes
GenomeDISCO|transition	yes
SGE|text	"-l h_vmem=3G"
slurm|text	"--mem 3G"
```
Note: all of the above parameters need to be specified in the parameters file.

Here are details about setting these parameters:

- `GenomeDISCO|subsampling` This allows subsampling the datasets to a specific desired sequencing depth. Possible values are: `lowest` (subsample to the depth of the sample with the lower sequencing depth from the pair being compared), `<samplename>` where <samplename> is the name of the sample that is used to determine the sequencing depth to subsample from. 

- `GenomeDISCO|tmin` The minimum number of steps of random walk to perform. Integer, > 0.

- `GenomeDISCO|tmax` The max number of steps of random walk to perform. Integer, > tmin.
 
- `GenomeDISCO|norm` The normalization to use on the data when running GenomeDISCO. Possible values include: `uniform` (no normalization), `sqrtvc`.

- `GenomeDISCO|scoresByStep` Whether to report the score at each t. By default (GenomeDISCO|scoresByStep no), only the final reproducibility score is returned.

- `GenomeDISCO|removeDiag` Whether to set the diagonal to entries in the contact map to 0. By default (GenomeDISCO|removeDiag yes), the diagonal entries are set to 0.

- `GenomeDISCO|transition` Whether to convert the normalized contact map to an appropriate transition matrix before running the random walks. By default (GenomeDISCO|transition yes) the normalized contact map is converted to a proper transition matrix, such that all rows sum to 1 exactly.

- `SGE|text` Text to append to the job submission for SGE. The default is "-l h_vmem=3G".

- `slurm|text` Text to append to the job submission for slurm. The default is "--mem 3G". 

Running GenomeDISCO step by step
============================================
GenomeDISCO consists of multiple steps, which are run in sequence by default. However, the user may decide to run the steps individually, which can be useful for instance when running GenomeDISCO with job submission engines that runs the comparisons in parallel as separate jobs.

**GenomeDISCO steps**

**preprocess**

Preprocesses all datasets provided in `--metadata_samples`.

Example command: 
```
genomedisco preprocess --metadata_samples examples/metadata.samples --bins examples/Bins.w50000.bed.gz --outdir examples/output --parameters_file examples/example_parameters.txt
```

**concordance**

Runs GenomeDISCO on all samples pairs provided in `--metadata_pairs`. 

Example command: 
```
genomedisco concordance --metadata_pairs examples/metadata.pairs --outdir examples/output 
```

**summary**

Summarizes scores across all comparisons.

Example command: 
```
genomedisco summary --metadata_samples examples/metadata.samples --metadata_pairs examples/metadata.pairs --bins examples/Bins.w50000.bed.gz --outdir examples/output 
```

**cleanup**

Clean up superfluous files, leaving only the scores.

Example command: 
```
genomedisco cleanup --outdir examples/output
```

Running GenomeDISCO with job submission engines
============

It is possible to run GenomeDISCO with job submission engines, specifically either SGE or slurm.
To do so, modify the parameters `SGE|text` or `slurm|text` respectively, to add any additional parameters to the job run.

Then, run the steps sequentially (that is, wait for all jobs of a given step to complete before launching the next step), while specifying `--running_mode` to either `sge` or `slurm`.

For instance, an example analysis workflow for SGE would be:
```
genomedisco preprocess --running_mode sge --metadata_samples examples/metadata.samples --bins examples/Bins.w50000.bed.gz --outdir examples/output --parameters_file examples/example_parameters.txt
genomedisco concordance --running_mode sge --metadata_pairs examples/metadata.pairs --outdir examples/output 
genomedisco summary --running_mode sge --metadata_samples examples/metadata.samples --metadata_pairs examples/metadata.pairs --bins examples/Bins.w50000.bed.gz --outdir examples/output 
genomedisco cleanup --running_mode sge --outdir examples/output
```

Similarly, for slurm, change sge to slurm for the `--running_mode`.

More questions?
====
Submit an issue for this repository.

This code was put together by Oana Ursu (oursu@stanford.edu).







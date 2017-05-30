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

Input files
---
In this code, we will refer to the contact map as a network, where nodes are genomic regions and edges are the potential contacts formed between the nodes.

To start, you will need the following:

**1. Contact maps to be compared**

2 contact map files for which you want to compute reproducibility. These should be in the format: chr1, node1, chr2, node2, value (tab separated). node1 and node2 represent 2 genomic regions and value represents the number of reads supporting the contact between node1 and node2.

An example line in this file:

`chr21 20000 chr21 150000 54`

**2. Nodes file**

A bed file that specifies the positions of the nodes, in the format: chr, start, end, name (tab separated). **Important** The name given to each node should correspond to what is used as node1 and node2 in the contact maps described above. For instance, in the example contact map we used above, the name is the start of the node, so an example line in the nodes file would be:

`chr21 20000 30000 20000`

However, some papers use the midpoint of the node, which would transform the example line above to:

`chr21 20000 30000 25000`

As long as the node name corresponds to what is used in the contact map, you are good to go.

Running GenomeDISCO
---

Say you want to compare 2 contact maps. For this example, we will use a subset of datasets from Rao et al., 2014. The data used below is in `genomedisco/genomedisco/examples`.

**1. Split files by chromosome**

```
splitByChromosome.sh -t hic -i examples/metadata.samples -n examples/Nodes.w40000.bed.gz -o examples/output
```

**2. Run GenomeDISCO**

```
genomedisco_GenomewideIntraChromosomal.sh -t hic -i examples/metadata.pairs -n examples/Nodes.w40000.bed.gz -o examples/output -b sqrtvc
```

Note: the output directory specified here should be the same as the one used in Step 1.

**3. Visualize results **

Create a beautiful html report describing the reproducibility analysis.

```
genomedisco_GenomewideIntraChromosomal_report.sh -t hic -i examples/metadata.pairs -n examples/Nodes.w40000.bed.gz -o examples/output -b sqrtvc
```

For the example we just ran, the html is here: http://htmlpreview.github.io/?http://github.com/kundajelab/genomedisco/blob/master/examples/output/results/sample1.vs.sample2/genomewide.sample1.vs.sample2.genomedisco.report.html

Running reproducibility analysis in batches
====
Want to run multiple comparisons, not just one? Here is how to modify the above commands.

**1. Split files by chromosome**

All you need to do here is add any additional samples involved in your comparisons.

```
HIC001 examples/HIC001.res40000.gz
HIC002 examples/HIC002.res40000.gz
HIC050 examples/HIC050.res40000.gz
```

Now, let's move on to split the data by chromosome.

```
cd genomedisco

#input files
nodes=$(pwd)/examples/Nodes.w40000.bed.gz
contactmap1=$(pwd)/examples/HIC001.res40000.gz
contactmap2=$(pwd)/examples/HIC002.res40000.gz
contactmap3=$(pwd)/examples/HIC050.res40000.gz ######## new addition

#create metadata for samples
metadata_samples=$(pwd)/examples/metadata.batch.samples
echo "sample1 ${contactmap1}" > ${metadata_samples}
echo "sample2 ${contactmap2}" >> ${metadata_samples}
echo "sample3 ${contactmap3}" >> ${metadata_samples} ######## new addition

#split input files by chromosome
outputdir=$(pwd)/examples/output_batch
scripts/splitByChromosome.sh -t hic -i ${metadata_samples} -n ${nodes} -o ${outputdir}
```

This will create a set of files in `$(pwd)/examples/output_batch`. These will be used in the next step.

**2. Run GenomeDISCO**

All you need to do here is add all comparisons to the metadata describing all pairs of datasets you want to compare.

```
cd genomedisco

#remember input data (this will be used to determine what chromosomes to compute)
nodes=$(pwd)/examples/Nodes.w40000.bed.gz

#create metadata for pairs to compare
metadata_pairs=$(pwd)/examples/metadata.batch.pairs
echo "sample1 sample2" > ${metadata_pairs}
echo "sample1 sample3" >> ${metadata_pairs} ######## new addition
echo "sample2 sample3" >> ${metadata_pairs} ######## new addition

#run reproducibility analysis
outputdir=$(pwd)/examples/output_batch
normalization=sqrtvc
scripts/genomedisco_GenomewideIntraChromosomal.sh -t hic -i ${metadata_pairs} -n ${nodes} -o ${outputdir}
-b ${normalization} 
```

**3. Generate genomewide report**

Nothing needs to be changed here, since this uses the metadata for pairs that we defined above.

```
cd genomedisco

#remember input data (this will be used to determine what chromosomes to compute)
nodes=$(pwd)/examples/Nodes.w40000.bed.gz
metadata_pairs=$(pwd)/examples/metadata.batch.pairs

#run reproducibility analysis
outputdir=$(pwd)/examples/output_batch
normalization=sqrtvc
scripts/genomedisco_GenomewideIntraChromosomal_report.sh -t hic -i ${metadata_pairs} -n ${nodes} -o ${outputdir} -b ${normalization}
```

Visualize the pretty reports as html files. Note that sample1 and sample2 are from the same cell type, while sample3 is a different cell type. This is reflected in the reproducibility scores. 
- sample1 vs sample2: http://htmlpreview.github.io/?http://github.com/kundajelab/genomedisco/blob/master/examples/output_batch/results/sample1.vs.sample2/genomewide.sample1.vs.sample2.genomedisco.report.html
- sample1 vs sample3: http://htmlpreview.github.io/?http://github.com/kundajelab/genomedisco/blob/master/examples/output_batch/results/sample1.vs.sample3/genomewide.sample1.vs.sample3.genomedisco.report.html
- sample2 vs sample3: http://htmlpreview.github.io/?http://github.com/kundajelab/genomedisco/blob/master/examples/output_batch/results/sample2.vs.sample3/genomewide.sample2.vs.sample3.genomedisco.report.html

More questions?
====
Contact Oana Ursu

oursu@stanford.edu

We are excited to hear from you and how we can improve GenomeDISCO!


# genomedisco
[![Build Status](https://travis-ci.org/kundajelab/genomedisco.svg?branch=master)](https://travis-ci.org/kundajelab/genomedisco)
[![license](https://img.shields.io/github/license/mashape/apistatus.svg?maxAge=2592000)](https://github.com/kundajelab/genomedisco/blob/master/LICENSE)

`genomedisco` is a package for comparing contact maps of 3D genome structures, obtained from experiments such as Hi-C, Capture-C, ChIA-PET, HiChip, etc. It uses graph diffusion to smooth contact maps, and then compares them, resulting in a reproducibility score that can be used for quality control of biological replicates.

Installation
===

1. Install [Anaconda](https://www.continuum.io/downloads). 
2. Install genomedisco with the following command:
```
conda install genomedisco -c kundajelab
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

**1. Split files by chromosome**

In addition to the input files, you need a metadata file with all the samples you are going to compare later on. The format is: contact map name, path to contact map (tab delimited). Here's an example:

```HIC001 examples/HIC001.res40000.gz
HIC002 examples/HIC002.res40000.gz
```

Now, let's move on to split the data by chromosome.

```
cd genomedisco

#input files
nodes=$(pwd)/examples/Nodes.w40000.bed.gz
contactmap1=$(pwd)/examples/HIC001.res40000.gz
contactmap2=$(pwd)/examples/HIC002.res40000.gz

#create metadata
metadata_samples=$(pwd)/examples/metadata.samples
echo "contactmapname1 ${contactmap1}" > ${metadata_samples}
echo "contactmapname2 ${contactmap2}" >> ${metadata_samples}

#split input files by chromosome
outputdir=$(pwd)/examples/output
scripts/splitByChromosome.sh -t hic -i ${metadata_samples} -n ${nodes} -o ${outputdir}
```

This will create a set of files in `${outputdir}/data`. These will be used in the next step.

**2. Run GenomeDISCO**



**3. Generate genomewide report**






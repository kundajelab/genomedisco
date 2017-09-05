#!/bin/sh
. /oak/stanford/groups/akundaje/oursu/code/new_genomedisco/genomedisco/scripts/bashrc.allMethods
${mypython} /oak/stanford/groups/akundaje/oursu/code/new_genomedisco/genomedisco/reproducibility_analysis/make_partition_from_bedfile.py --nodes /oak/stanford/groups/akundaje/oursu/code/new_genomedisco/genomedisco/examples/Nodes.w40000.bed.gz --partition /oak/stanford/groups/akundaje/oursu/code/new_genomedisco/genomedisco/examples/output/data/forQuASAR/nodes.partition --subset_chromosomes NA --resolution 40000

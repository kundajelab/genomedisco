#!/bin/sh
. /oak/stanford/groups/akundaje/oursu/code/new_genomedisco/genomedisco/scripts/bashrc.allMethods
zcat -f /oak/stanford/groups/akundaje/oursu/code/new_genomedisco/genomedisco/examples/Nodes.w40000.bed.gz | sort -k1,1 -k2,2n | awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\tincluded"}' | sed 's/chrchr/chr/g' | awk -v chromosome=chr21 '{if ($1==chromosome) print $0}' | gzip > /oak/stanford/groups/akundaje/oursu/code/new_genomedisco/genomedisco/examples/output/data/nodes/nodes.chr21.gz

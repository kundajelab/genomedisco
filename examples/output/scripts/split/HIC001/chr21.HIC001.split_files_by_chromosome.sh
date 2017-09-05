#!/bin/sh
. /oak/stanford/groups/akundaje/oursu/code/new_genomedisco/genomedisco/scripts/bashrc.allMethods
mkdir -p /oak/stanford/groups/akundaje/oursu/code/new_genomedisco/genomedisco/examples/output/data/edges/HIC001
zcat -f /oak/stanford/groups/akundaje/oursu/code/new_genomedisco/genomedisco/examples/HIC001.res40000.gz | awk '{print "chr"$1"\t"$2"\tchr"$3"\t"$4"\t"$5}' | sed 's/chrchr/chr/g' | awk -v chromosome=chr21 '{if ($1==chromosome && $3==chromosome) print $2"\t"$4"\t"$5}' | gzip > /oak/stanford/groups/akundaje/oursu/code/new_genomedisco/genomedisco/examples/output/data/edges/HIC001/HIC001.chr21.gz

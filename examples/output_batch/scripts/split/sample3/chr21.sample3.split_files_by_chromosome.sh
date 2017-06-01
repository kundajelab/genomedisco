#!/bin/sh
source /srv/gsfs0/projects/snyder/oursu/software/git/public_genomedisco/genomedisco/scripts/bashrc_genomedisco
mkdir -p /srv/gsfs0/projects/snyder/oursu/software/git/public_genomedisco/genomedisco/examples/output_batch/data/edges/sample3
zcat -f /srv/gsfs0/projects/snyder/oursu/software/git/public_genomedisco/genomedisco/examples/HIC050.res40000.gz | awk '{print "chr"$1"\t"$2"\tchr"$3"\t"$4"\t"$5}' | sed 's/chrchr/chr/g' | awk -v chromosome=chr21 '{if ($1==chromosome && $3==chromosome) print $2"\t"$4"\t"$5}' | gzip > /srv/gsfs0/projects/snyder/oursu/software/git/public_genomedisco/genomedisco/examples/output_batch/data/edges/sample3/sample3.chr21.gz

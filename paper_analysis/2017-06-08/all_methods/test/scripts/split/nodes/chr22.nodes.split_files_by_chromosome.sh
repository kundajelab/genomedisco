#!/bin/sh
source /srv/gsfs0/projects/snyder/oursu/software/git/public_genomedisco/genomedisco/scripts/bashrc_genomedisco
zcat -f /srv/gsfs0/projects/snyder/oursu/software/git/public_genomedisco/genomedisco/paper_analysis/2017-02-05/all_methods/test/prefix.bins.gz | awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\tincluded"}' | sed 's/chrchr/chr/g' | awk -v chromosome=chr22 '{if ($1==chromosome) print $0}' | gzip > /srv/gsfs0/projects/snyder/oursu/software/git/public_genomedisco/genomedisco/paper_analysis/2017-02-05/all_methods/test/data/nodes/nodes.chr22.gz

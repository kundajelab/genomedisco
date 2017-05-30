#!/bin/sh
source scripts/bashrc_genomedisco
zcat -f /srv/gsfs0/projects/snyder/oursu/software/git/public_genomedisco/genomedisco/examples/HIC002.res40000.gz | awk '{print "chr"$1"\t"$2"\tchr"$3"\t"$4"\t"$5}' | sed 's/chrchr/chr/g' | awk -v chromosome=chr21 '{if ($1==chromosome && $3==chromosome) print $2"\t"$4"\t"$5}' | gzip > /srv/gsfs0/projects/snyder/oursu/software/git/public_genomedisco/genomedisco/examples/output_batch/data/edges/sample2/sample2.chr21.gz
echo DONE
zcat -f /srv/gsfs0/projects/snyder/oursu/software/git/public_genomedisco/genomedisco/examples/HIC002.res40000.gz | awk '{print "chr"$1"\t"$2"\tchr"$3"\t"$4"\t"$5}' | sed 's/chrchr/chr/g' | awk -v chromosome=chr22 '{if ($1==chromosome && $3==chromosome) print $2"\t"$4"\t"$5}' | gzip > /srv/gsfs0/projects/snyder/oursu/software/git/public_genomedisco/genomedisco/examples/output_batch/data/edges/sample2/sample2.chr22.gz
echo DONE

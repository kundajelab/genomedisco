
import argparse
import copy
import re
import os
from time import gmtime, strftime
import gzip
import numpy as np

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--scores_by_chromosome',default='/ifs/scratch/oursu/paper_2017-12-20/results/rao/res50000.final/compiled_scores/GenomeDISCO.scores.multiple_t.by_chromosome.txt.gz')
    parser.add_argument('--out',default='/ifs/scratch/oursu/paper_2017-12-20/results/rao/res50000.final/compiled_scores/GenomeDISCO.scores.multiple_t.genomewide.txt.gz')
    args = parser.parse_args()

    score_dict={}
    for line in gzip.open(args.scores_by_chromosome):
        items=line.strip().split('\t')
        chromo,m1,m2=items[0],items[1],items[2]
        print chromo
        comparison=m1+'\t'+m2
        if comparison not in score_dict:
            score_dict[comparison]={}

        score_cols=range(3,len(items))
        for score_col in score_cols:
            if score_col not in score_dict[comparison]:
                score_dict[comparison][score_col]=[]
            score=float(items[score_col])
            score_dict[comparison][score_col].append(score)
    
    out=gzip.open(args.out,'w')
    for comparison in score_dict:
        scores=[]
        for score_col_idx in range(len(score_cols)):
            score_col=score_cols[score_col_idx]
            genomewide=np.mean(np.array(score_dict[comparison][score_col]))
            scores.append(str(genomewide))
        out.write('genomewide'+'\t'+comparison+'\t'+'\t'.join(scores)+'\n')
    out.close()
    print args.out

main()

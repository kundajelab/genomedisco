
import argparse
import numpy as np

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--metadata_pairs',default='/srv/gsfs0/projects/kundaje/users/oursu/3d/encode_highres/AllChrAnon/processed/metadata/metadata.res40000.pairs')
    parser.add_argument('--scoring_file',default='/ifs/scratch/oursu/paper_2017-12-20/encode_highres.final/res40000/results/compiled_scores.txt')
    parser.add_argument('--out',default='test')
    parser.add_argument('--chromo_order',default='chr1,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr2,chr20,chr21,chr22,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chrX')

    args = parser.parse_args()

    chromos=args.chromo_order.split(',')
    out=open(args.out,'w')
    out.write('#m1'+'\t'+'m2'+'\t'+'\t'.join(chromos)+'\t'+'genomewide'+'\n')

    scores={}
    for line in open(args.scoring_file,'r').readlines():
        items=line.strip().split('\t')
        m1,m2,chromo,score=items[0],items[1],items[2],float(items[3])
        if m1 not in scores:
            scores[m1]={}
        if m2 not in scores:
            scores[m2]={}
        if m1 not in scores[m2]:
            scores[m2][m1]={}
        if m2 not in scores[m1]:
            scores[m1][m2]={}
        if chromo not in scores[m1][m2]:
            scores[m1][m2][chromo]=score
        if chromo not in scores[m2][m1]:
            scores[m2][m1][chromo]=score

    seen=set()
    for line in open(args.metadata_pairs,'r').readlines():
        items=line.strip().split('\t')
        m1,m2=items[0],items[1]
        if m1+'vs'+m2 in seen:
            continue
        seen.add(m1+'vs'+m2)
        stuff=[m1,m2]
        scorelist=[]
        for chromo in chromos:
            cur_score=scores[m1][m2][chromo]
            scorelist.append(cur_score)
            stuff.append(str(cur_score))
        stuff.append(str(np.mean(np.array(scorelist))))
        out.write('\t'.join(stuff)+'\n')
        
    

main()

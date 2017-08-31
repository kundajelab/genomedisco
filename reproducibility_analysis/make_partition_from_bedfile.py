import argparse
import sys
import os

import hifive
import h5py
import numpy
import gzip
import re
import subprocess as subp

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--nodes')
    parser.add_argument('--partition')
    parser.add_argument('--resolution',type=int)
    parser.add_argument('--re',action='store_true')
    parser.add_argument('--subset_chromosomes',default='NA')
    args = parser.parse_args()

    nodes_modified_file=args.partition+'.tmp'
    subp.check_output(['bash','-c','zcat -f '+args.nodes+' | sort -k1,1 -k2,2n | gzip > '+nodes_modified_file+'.sorted'])
    nodes_modified=open(nodes_modified_file,'w')
    for line in gzip.open(nodes_modified_file+'.sorted','r'):
        items=line.strip().split('\t')
        chromo,start,end,name=items[0],items[1],items[2],items[3]
        if args.subset_chromosomes!='NA':
            if chromo not in args.subset_chromosomes.split(','):
                continue
        nodes_modified.write(re.sub('chr','',chromo)+'\t'+start+'\t'+end+'\t'+name+'\n')
    nodes_modified.close()

    if args.re:
        #restriction fragments
        myfends=hifive.Fend(args.partition, mode='w')
        myfends.load_fends(nodes_modified_file, format='bed') 
    else:
        #uniform bins
        myfends=hifive.Fend(args.partition, mode='w',binned=int(args.resolution)) 
        myfends.load_bins(nodes_modified_file, format='bed') 
    myfends.save() 

    os.remove(nodes_modified_file)
    os.remove(nodes_modified_file+'.sorted')

main()

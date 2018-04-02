
import argparse
import copy
import re
import os
from time import gmtime, strftime

from genomedisco import data_operations, processing, visualization
from genomedisco.comparison_types.disco_random_walks import DiscoRandomWalks
from genomedisco.comparison_types.disco_random_walks import to_transition

def main():
    parser = argparse.ArgumentParser(description='Compute RW transformation of 3D data')
    parser.add_argument('--datatype',default='hic')
    parser.add_argument('--m',type=str)
    parser.add_argument('--matrix_format',type=str,default='n1n2val',help='c1n1c2n2val')
    parser.add_argument('--node_file',type=str)
    parser.add_argument('--remove_diagonal',action='store_true')
    parser.add_argument('--mname',type=str)
    parser.add_argument('--outdir',type=str,default='OUT')
    parser.add_argument('--outpref',type=str,default='outpref')
    parser.add_argument('--norm',type=str,default='uniform')
    parser.add_argument('--method',type=str,default='RandomWalks')
    parser.add_argument('--tmin',type=int,default=1)
    parser.add_argument('--tmax',type=int,default=3)
    parser.add_argument('--transition',action='store_true')
    parser.add_argument('--blacklist',default='NA')
    args = parser.parse_args()

    os.system('mkdir -p '+args.outdir)
    nodes,nodes_idx,blacklist_nodes=processing.read_nodes_from_bed(args.node_file,args.blacklist)

    m=processing.construct_csr_matrix_from_data_and_nodes(args.m,nodes,blacklist_nodes,args.remove_diagonal)

    m_norm=data_operations.process_matrix(m,args.norm)

    mup=m_norm
    mdown=mup.transpose()
    mdown.setdiag(0)
    m_full=mup+mdown

    if args.transition:
        m_full=to_transition(m_full)

    outname=args.outdir+'/'+args.outpref
    for t in range(args.tmin,(args.tmax+1)):
        if t==1:
            rw=copy.deepcopy(m_full)
        else:
            rw=rw.dot(m_full)
        processing.write_matrix_from_csr_and_nodes(rw,nodes_idx,outname+'.rw_t'+str(t)+'.gz')


if __name__=="__main__":
    main()

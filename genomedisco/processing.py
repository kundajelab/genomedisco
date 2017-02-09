import numpy as np
import gzip
from scipy.sparse import csr_matrix
import logging

#===== MATRIX IO
#from http://stackoverflow.com/questions/8955448/save-load-scipy-sparse-csr-matrix-in-portable-data-format
def save_sparse_csr(filename,array):
    np.savez(filename,data = array.data ,indices=array.indices,
             indptr =array.indptr, shape=array.shape )

def load_sparse_csr(filename):
    loader = np.load(filename)
    return csr_matrix((  loader['data'], loader['indices'], loader['indptr']),
                         shape = loader['shape'])

#===== Reading in data
#TODO: tell people that nodes should come in the order in which they go in the matrix
def read_nodes_from_bed(bedfile):
    logging.info("| processing: Loading genomic regions from "+bedfile)

    nodes={}
    nodes_idx={}
    node_c=0
    for line in gzip.open(bedfile,'r'):
        items=line.strip().split('\t')
        chromo=items[0]
        start=items[1]
        end=items[2]
        node=items[3]
        if len(items)>4:
            include=items[4]
        
        if node in nodes.keys():
            logging.error("Error: Genomic region appears multiple times in your file. One such example is "+node+". Please make sure all genomic regions are unique and re-run")
            sys.exit()
        if node not in nodes.keys():
            nodes[node]={}
            nodes[node]['idx']=node_c
            nodes[node]['chr']=chromo
            nodes[node]['start']=start
            nodes[node]['end']=end
            if len(items)>4:
                nodes[node]['include']=include
            nodes_idx[node_c]=node 
            node_c+=1
    return nodes,nodes_idx

def construct_csr_matrix_from_data_and_nodes(f,nodes,remove_diag=True):
    logging.info("| processing: Loading interaction data from "+f)

    total_nodes=len(nodes.keys())
    mdata=np.loadtxt(f)
    
    i=map(lambda x:nodes[str(int(x))]['idx'], mdata[:,0])
    j=map(lambda x:nodes[str(int(x))]['idx'], mdata[:,1])
    
    #flag cases where the same i,j pair is repeated in the file                                                                                                            
    #- convert i,j,value to min(i,j),max(i,j),value                                                                                                                        
    ij=np.array([i,j])
    mini=ij.min(axis=0)
    maxi=ij.max(axis=0)
    mini_maxi_ij=np.array([mini,maxi]).T
    #- have a list of rows, where each row is a tuple of 2 nodes
    rows=[tuple(row) for row in mini_maxi_ij]
    #- if the original set of rows is larger than the unique set of rows, flag an error
    if len(rows)>len(set(rows)):
        logging.warning("| processing: =============== Warning: Your file contains duplicate interactions! Please ensure that each interaction is listed once, then re-run. In the meantime, we will run this analysis using the sum of all counts encountered per interaction")
    
    csr_m=csr_matrix( (mdata[:,2],(mini_maxi_ij[:,0],mini_maxi_ij[:,1])), shape=(total_nodes,total_nodes),dtype=float )
    if remove_diag:
        csr_m.setdiag(0)
    
    return csr_m


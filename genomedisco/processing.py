import numpy as np
import gzip
from scipy.sparse import csr_matrix
from scipy.sparse import coo_matrix
from time import gmtime, strftime

#===== MATRIX IO
#from http://stackoverflow.com/questions/8955448/save-load-scipy-sparse-csr-matrix-in-portable-data-format
def save_sparse_csr(filename,array):
    np.savez(filename,data = array.data ,indices=array.indices,
             indptr =array.indptr, shape=array.shape )

def load_sparse_csr(filename):
    loader = np.load(filename)
    return csr_matrix((  loader['data'], loader['indices'], loader['indptr']),
                         shape = loader['shape'])

def read_nodes_from_bed(bedfile,blacklistfile='NA'):
    
    blacklist={}
    if blacklistfile!='NA':
        for line in gzip.open(blacklistfile):
            items=line.strip().split('\t')
            chromo,start,end=items[0],int(items[1]),int(items[2])
            if chromo not in blacklist:
                blacklist[chromo]=[]
            blacklist[chromo].append((start,end))
    
    print "GenomeDISCO | "+strftime("%c")+" | processing: Loading genomic regions from "+bedfile

    nodes={}
    nodes_idx={}
    node_c=0
    blacklisted_nodes=[]
    for line in gzip.open(bedfile,'r'):
        items=line.strip().split('\t')
        chromo=items[0]
        start=int(items[1])
        end=int(items[2])
                
        node=items[3]
        if len(items)>4:
            include=items[4]
        
        if node in nodes.keys():
            print "GenomeDISCO | "+strftime("%c")+" | Error: Genomic region appears multiple times in your file. One such example is "+node+". Please make sure all genomic regions are unique and re-run"
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
            
            if chromo in blacklist:
                for blacklist_item in blacklist[chromo]:
                    if (start<=blacklist_item[0] and end>=blacklist_item[0]) or (start<=blacklist_item[1] and end>=blacklist_item[1]) or (start>=blacklist_item[0] and end<=blacklist_item[1]):
                        blacklisted_nodes.append(node_c)
                
            node_c+=1
            
    return nodes,nodes_idx,blacklisted_nodes

def filter_nodes(m,to_remove):
    
    if len(to_remove)==0:
        return m
    
    nonzeros=m.nonzero()
    num_elts=len(nonzeros[0])
    
    r_idx=[i for i, x in enumerate(nonzeros[0]) if x not in to_remove]
    c_idx=[i for i, x in enumerate(nonzeros[1]) if x not in to_remove]
    keep=list(set(r_idx).union(set(c_idx)))
    
    coo_mat=m.tocoo()
        
    return csr_matrix((coo_mat.data[keep],(coo_mat.row[keep],coo_mat.col[keep])),shape=m.get_shape(),dtype=float) 
    

def construct_csr_matrix_from_data_and_nodes(f,nodes,blacklisted_nodes=[],remove_diag=True):
    print "GenomeDISCO | "+strftime("%c")+" | processing: Loading interaction data from "+f

    total_nodes=len(nodes.keys())
    i=[]
    j=[]
    v=[]

    #print strftime("%c")
    c=0
    for line in gzip.open(f):
        items=line.strip().split('\t')
        n1,n2,val=nodes[items[0]]['idx'],nodes[items[1]]['idx'],float(items[2])
        mini=min(n1,n2)
        maxi=max(n1,n2)
        i.append(mini)
        j.append(maxi)
        v.append(val)
        c+=1

    csr_m=csr_matrix( (v,(i,j)), shape=(total_nodes,total_nodes),dtype=float)
    if remove_diag:
        csr_m.setdiag(0)
    return filter_nodes(csr_m,blacklisted_nodes)

def old_construct_csr_matrix_from_data_and_nodes(f,nodes,blacklisted_nodes,remove_diag=True):
    print "GenomeDISCO | "+strftime("%c")+" | processing: Loading interaction data from "+f

    total_nodes=len(nodes.keys())
    mdata=np.loadtxt(f)
    
    dist_threshold=2000000

    #keep=(abs(mdata[:,0]-mdata[:,1])<=dist_threshold)
    #mdata=mdata[keep,:]
    
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
        print "=============== Warning: Your file contains duplicate interactions! Please ensure that each interaction is listed once, then re-run. In the meantime, we will run this analysis using the sum of all counts encountered per interaction"
    
    csr_m=csr_matrix( (mdata[:,2],(mini_maxi_ij[:,0],mini_maxi_ij[:,1])), shape=(total_nodes,total_nodes),dtype=float )
    if remove_diag:
        csr_m.setdiag(0)
    
    return filter_nodes(csr_m,blacklisted_nodes)




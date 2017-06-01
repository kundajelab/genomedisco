import sys
import copy
import random
import numpy as np
from scipy.sparse import csr_matrix
import math
from time import gmtime, strftime

#todo: add bait vs not bait information
def get_distance_dep_using_nodes_capturec(m,nodes,nodes_idx,approximation=10000):
    assert m.shape[0]==m.shape[1]
    dcounts={}
    pcounts={}
    total_reads=0.0
    marray=copy.deepcopy(m).toarray()
    entries_per_key={}
    
    include=set()
    for node in nodes:
        if nodes[node]['include']=='included':
            include.add(nodes[node]['idx'])
    for i in include:
        for j in range(m.shape[0]):
            if j in include:
                continue
            v=m[i,j]
            istart=int(nodes[nodes_idx[i]]['start'])
            jstart=int(nodes[nodes_idx[j]]['start'])
            di=math.ceil(1.0*abs(istart-jstart)*1.0/approximation)
            if di not in dcounts:
                dcounts[di]=0.0
                pcounts[di]=0.0
                entries_per_key[di]=0.0
            dcounts[di]+=v
            entries_per_key[di]+=1
            total_reads+=m[i,j]
    for di in dcounts:
        pcounts[di]=1.0*dcounts[di]/(total_reads*entries_per_key[di])
    return pcounts

def get_distance_dep(m):
    assert m.shape[0]==m.shape[1]
    dcounts={}
    pcounts={}
    #get all distances we'll want to take into account
    for di in range(m.shape[0]):
        dcounts[di]=0
        pcounts[di]=0
    nonzeros=m.nonzero()
    num_elts=len(nonzeros[0])
    total_reads=0
    for elt in range(num_elts):
        i=nonzeros[0][elt]
        j=nonzeros[1][elt]
        di=abs(i-j)
        dcounts[di]+=m[i,j]
        total_reads+=m[i,j]
    for di in range(m.shape[0]):
        pcounts[di]=1.0*dcounts[di]/((m.shape[0]-di)*total_reads)
    return pcounts

def sqrtvc(m):
    mup=m
    mdown=mup.transpose()
    mdown.setdiag(0)
    mtogether=mup+mdown
    sums_sq=np.sqrt(mtogether.sum(axis=1)) 
    nonzeros=m.nonzero()
    num_elts=len(nonzeros[0])
    rows=[]
    cols=[]
    m_norm_data=[]
    for elt in range(num_elts):
        i=nonzeros[0][elt]
        j=nonzeros[1][elt]
        rows.append(i)
        cols.append(j)
        if sums_sq[i,0]>0 and sums_sq[j,0]>0:
            m_norm_data.append(float(m[i,j])/(float(sums_sq[i,0])*float(sums_sq[j,0])))
        else:
            m_norm_data.append(0)
    return csr_matrix((m_norm_data,(rows,cols)),shape=m.get_shape(),dtype=float)

def coverage_norm(m):
    mup=m
    mdown=mup.transpose()
    mdown.setdiag(0)
    mtogether=mup+mdown
    sums=mtogether.sum(axis=1)
    nonzeros=m.nonzero()
    num_elts=len(nonzeros[0])
    rows=[]
    cols=[]
    m_norm_data=[]
    for elt in range(num_elts):
        i=nonzeros[0][elt]
        j=nonzeros[1][elt]
        rows.append(i)
        cols.append(j)
        if sums[i,0]>0 and sums[j,0]>0:
            m_norm_data.append(float(m[i,j])/(float(sums[i,0])*float(sums[j,0])))
        else:
            m_norm_data.append(0)
    return csr_matrix((m_norm_data,(rows,cols)),shape=m.get_shape(),dtype=float)

#assumes matrix is upper triangular
def matrix_2_coverageVector(m):
    mup=m
    mdown=mup.transpose()
    mdown.setdiag(0)
    mtogether=mup+mdown
    sums=mtogether.sum(axis=1)
    return sums

def array_2_coverageVector(m):
    assert np.allclose(m, np.triu(m))
    m_sym=m+m.T-m.diagonal()
    return m_sym.sum(axis=0)

def uniform_processing(m):
    return m

def process_matrix(m,matrix_processing):
    if matrix_processing=='uniform':
        return uniform_processing(m)
    if matrix_processing=='coverage_norm':
        return coverage_norm(m)
    if matrix_processing=='sqrtvc':
        return sqrtvc(m)

def subsample_to_depth(m,seq_depth):
    if type(m) is csr_matrix:
        return subsample_to_depth_csr_upperTri(m,seq_depth)
    if type(m) is np.ndarray:
        return subsample_to_depth_array_upperTri(m,seq_depth)

def subsample_to_depth_array_upperTri(m,seq_depth):
    m=np.triu(m)
    subsampled_data=np.zeros(m.shape)
    depthm=m.sum()
    assert seq_depth<=depthm
    subsampling_prob=seq_depth/depthm
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            if j<=i:
                continue
            n=m[i,j]
            subsampled_data[i,j]=np.random.binomial(n,subsampling_prob,1)[0]
    return subsampled_data

def subsample_to_depth_csr_upperTri(m,seq_depth):
    depthm=m.sum()
    assert seq_depth<=depthm
    subsampling_prob=seq_depth/depthm
    nonzeros=m.nonzero()
    num_elts=len(nonzeros[0])
    rows=[]
    cols=[]
    m_subsampled_data=[]
    for elt in range(num_elts):
        i=nonzeros[0][elt]
        j=nonzeros[1][elt]
        rows.append(i)
        cols.append(j)
        n=m[i,j]
        m_subsampled_data.append(np.random.binomial(n,subsampling_prob,1)[0])
    return csr_matrix((m_subsampled_data,(rows,cols)),shape=m.get_shape(),dtype=float)



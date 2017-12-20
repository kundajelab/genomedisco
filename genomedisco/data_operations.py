import sys
import copy
import random
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse import coo_matrix
import scipy as scipy
import math
from time import gmtime, strftime
import cProfile
import timeit
import scipy.sparse as sps

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
    D_sq = sps.spdiags(1.0/sums_sq.flatten(), [0], mtogether.get_shape()[0], mtogether.get_shape()[1], format='csr')
    return sps.triu(D_sq.dot(mtogether.dot(D_sq)))

def hichip_add_diagonal(m):
    mup=m
    mdown=mup.transpose()
    mdown.setdiag(0)
    mtogether=mup+mdown
    sums=mtogether.sum(axis=1)
    max_sum=np.max(sums)
    to_add=1.0*max_sum-1.0*sums
    to_add_values=[]
    for i in range(m.shape[0]):
        to_add_values.append(to_add[i,0])
    mtogether.setdiag(np.array(to_add_values))
    D = sps.spdiags(1.0/sums.flatten(), [0], mtogether.get_shape()[0], mtogether.get_shape()[1], format='csr')
    return sps.triu(D.dot(mtogether))

def coverage_norm(m):
    mup=m
    mdown=mup.transpose()
    mdown.setdiag(0)
    mtogether=mup+mdown
    sums=mtogether.sum(axis=1)
    D = sps.spdiags(1.0/sums.flatten(), [0], mtogether.get_shape()[0], mtogether.get_shape()[1], format='csr')
    return sps.triu(D.dot(mtogether.dot(D)))

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
    if matrix_processing=='fill_diagonal':
        return hichip_add_diagonal(m)

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

    vals=m.data
    num_elts=len(vals)
    m_subsampled_data=[]#np.random.binomial(value,subsampling_prob)
    elt=0
    while elt<num_elts:
        m_subsampled_data.append(np.random.binomial(vals[elt],subsampling_prob,1)[0])
        elt+=1
    return csr_matrix((m_subsampled_data, m.indices, m.indptr), dtype=int,shape=m.shape)
    

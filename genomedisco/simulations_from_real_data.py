import argparse
import copy
import re
import os
import gzip
from time import gmtime, strftime
import numpy as np
from scipy.sparse import csr_matrix

from genomedisco import data_operations, processing, visualization

def main():
    parser = argparse.ArgumentParser(description='Simulate Hi-C data based on real datasets.')
    parser.add_argument('--matrices',default='/ifs/scratch/oursu/3d/paper/2017-06-08/LA/reproducibility/res40000/data/edges/HIC001/HIC001.chr21.gz')
    parser.add_argument('--matrix_names',default='HIC001')
    parser.add_argument('--nodes',default='/ifs/scratch/oursu/3d/paper/2017-06-08/LA/reproducibility/res40000/data/nodes/nodes.chr21.gz')
    parser.add_argument('--distDepData',default='/ifs/scratch/oursu/3d/paper/2017-06-08/LA/reproducibility/res40000/data/edges/HIC001/HIC001.chr21.gz,/ifs/scratch/oursu/3d/paper/2017-06-08/LA/reproducibility/res40000/data/edges/HIC002/HIC002.chr21.gz')
    parser.add_argument('--edgenoise',default='0.0')
    parser.add_argument('--nodenoise',default='0.0')
    parser.add_argument('--boundarynoise',default='0')
    parser.add_argument('--depth',type=int,default=1000000)
    parser.add_argument('--outdir',default='/ifs/scratch/oursu/test/testmat')
    parser.add_argument('--resolution',type=int,default=40000)
    parser.add_argument('--mini',type=int,default=-1)
    parser.add_argument('--maxi',type=int,default=-1)
    args = parser.parse_args()

    #setup nodes
    nodes,nodes_idx,blacklisted_nodes=processing.read_nodes_from_bed(args.nodes)

    #set mini and maxi coordinates to focus on when simulating
    if args.mini<=-1:
        args.mini=0
    if args.maxi<=-1:
        args.maxi=len(nodes.keys())

    #now go through each of the matrices, and simulate from them
    matrices=args.matrices.split(',')
    matrix_names=args.matrix_names.split(',')
    for m_idx in range(len(matrices)):
        #read in the matrix
        mname=matrix_names[m_idx]
        mfile=matrices[m_idx]
        my_matrix_orig=read_in_data(mfile,nodes)
        for edgenoise in args.edgenoise.split(','):
            for nodenoise in args.nodenoise.split(','):
                for boundarynoise in args.boundarynoise.split(','):
                    my_matrix=shift_dataset(my_matrix_orig,int(boundarynoise))
                    ddfiles=args.distDepData.split(',')
                    for ddfile_idx in range(len(ddfiles)):
                        ddfile=ddfiles[ddfile_idx]
                        dd=read_in_data(ddfile,nodes)
                        prob_matrix=get_probability_matrix(my_matrix,dd,float(edgenoise),float(nodenoise),args.mini,args.maxi,np.random.RandomState(hash('probability')%10000))
                        ablist=['a','b']
                        for ab_idx in range(len(ablist)):
                            ab=ablist[ab_idx]
                            if ab=='a':
                                s=hash(mname+ddfile)%10000
                            if ab=='b':
                                s=hash(mname+ddfile)%10000+101
                            
                            intro=args.outdir+'/Depth_'+str(args.depth)+'.'+mname
                            sampled_matrix=sample_interactions(copy.deepcopy(prob_matrix),args.depth,np.random.RandomState(s))
                            ftowrite=intro+'.EN_'+str(edgenoise)+'.NN_'+str(nodenoise)+'.BN_'+str(boundarynoise)+'.'+ab+'.dd_'+str(ddfile_idx)+'.gz'
                            print ftowrite
                            write_matrix(sampled_matrix,ftowrite,args)
                            
def sample_interactions(prob_matrix1,depth,pet_random):
    new_m=np.zeros(prob_matrix1.shape)
    for i in range(new_m.shape[0]):
        for j in range(i,new_m.shape[0]):
            p=prob_matrix1[i,j]
            reads=pet_random.binomial(depth, p, size=1)[0]
            new_m[i,j]=reads
            new_m[j,i]=reads
    return new_m

def get_probability_matrix(my_matrix,ddmat,edge_noise,node_noise,mini,maxi,pet_random):#,maxdist=2000):

    #0 out diagonal of our matrix #==========
    for i in range(my_matrix.shape[0]):
        my_matrix[i,i]=0.0

    #convert matrix to probabilities
    total_probs=np.triu(copy.deepcopy(my_matrix)).sum()
    mat=copy.deepcopy(my_matrix)/total_probs #this is the probability matrix
    mat_dd=data_operations.get_distance_dep(csr_matrix(mat))
    #0 out diagonal of the ddmat #================
    for i in range(ddmat.shape[0]):
        ddmat[i,i]=0.0
    desired_dd=data_operations.get_distance_dep(csr_matrix(ddmat))
    
    #rescale the values to obey the distance curve given
    new_mat=np.zeros(mat.shape)
    mat_ddsums={}
    desired_ddsums={}
    mat_total=0.0
    desired_total=0.0
    for i in range(new_mat.shape[0]-1):
        mat_ddsums[i]=np.diagonal(mat,i).sum()
        desired_ddsums[i]=np.diagonal(ddmat,i).sum()
        mat_total+=mat_ddsums[i]
        desired_total+=desired_ddsums[i]
    for i in range(new_mat.shape[0]-1):
        for j in range(i,new_mat.shape[0]):#,i+maxdist+1)):
            d=abs(i-j)   
            if d not in mat_ddsums:
                continue
            if desired_total*mat_ddsums[d]==0.0:
                val=0.0
            else:
                val=1.0*mat[i,j]*mat_total*desired_ddsums[d]/(desired_total*mat_ddsums[d])
            new_mat[i,j]=val
            new_mat[j,i]=val

        
    #0 out things that are not within mini<->maxi
    for i in range(new_mat.shape[0]):
        if i>=mini and i<=maxi:
            continue
        new_mat[i,:]=0.0
        new_mat[:,i]=0.0        

    #rescale the new matrix to be a probability matrix
    total_probs=np.triu(copy.deepcopy(new_mat)).sum()
    new_mat=np.triu(copy.deepcopy(new_mat))/total_probs
    

    #edge noise
    #pet_random = np.random.RandomState()
    new_mat2=np.zeros(new_mat.shape)
    new_mat2=new_mat
    
    for i in range(new_mat.shape[0]):
        for j in range(min(i,new_mat.shape[0]),new_mat.shape[0]): 
            pij=new_mat[i,j]
            pij_final=pij
            if edge_noise!=0.0:
                boink = float(pet_random.binomial(1, edge_noise, size=1)[0]) #(pet_random.rand(1)[0] <= edge_noise)
                if (boink>0.0):
                    pij_final=0.0
            new_mat2[i,j]=pij_final
    
    #node noise
    for i in range(new_mat2.shape[0]):
        remove_node=float(pet_random.binomial(1, node_noise, size=1)[0])
        if remove_node>0.0:
            new_mat2[i,:]=0.0
            new_mat2[:,i]=0.0
    total_probs=np.triu(copy.deepcopy(new_mat2)).sum()
    return np.triu(new_mat2)/total_probs
    
def read_in_data(mname_full,nodes):
    mat=processing.construct_csr_matrix_from_data_and_nodes(mname_full,nodes,[],True).toarray()
    mat=mat + mat.T
    return mat

def write_matrix(sampled_matrix,fname,args,chromo='chr21'):
    mini=args.mini
    maxi=args.maxi
    f=gzip.open(fname,'w')
    for i in range(mini,min(maxi,sampled_matrix.shape[0])):
        for j in range(min(i+1,sampled_matrix.shape[0]),min(maxi,sampled_matrix.shape[0])):
            v=sampled_matrix[i,j]
            if v>0.0:
                n1=i*args.resolution
                n2=j*args.resolution
                f.write(chromo+'\t'+str(n1)+'\t'+chromo+'\t'+str(n2)+'\t'+str(v)+'\n')
    f.close()

def shift_dataset(m,boundarynoise):
    if boundarynoise==0:
        return m
    nonzero_rows=np.where(m.any(axis=1))[0]
    small_m=copy.deepcopy(m)
    small_m=small_m[nonzero_rows,:]
    small_m=small_m[:,nonzero_rows]
    print small_m
    print 'roll'
    small_m=np.roll(small_m,boundarynoise,axis=0)
    print small_m
    print 'roll2'
    small_m=np.roll(small_m,boundarynoise,axis=1)
    print small_m
    outm=np.zeros(m.shape)
    for i_idx in range(len(nonzero_rows)):
        i=nonzero_rows[i_idx]
        for j_idx in range(i_idx,len(nonzero_rows)):
            j=nonzero_rows[j_idx]
            outm[i,j]=small_m[i_idx,j_idx]
            outm[j,i]=outm[i,j]
    return outm

main()

        


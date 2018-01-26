import sys
import matplotlib
matplotlib.use('Agg')
import gzip
import numpy as np
import matplotlib.pyplot as plt
import os
import re
import copy
from scipy.stats.mstats import mquantiles
from scipy.spatial.distance import euclidean
from sklearn import metrics
from pylab import rcParams
from time import gmtime, strftime
from scipy.sparse import csr_matrix
from scipy import sparse
import h5py
import psutil
import scipy.sparse as sps

def to_transition(mtogether):
    sums=mtogether.sum(axis=1)
    #make the ones that are 0, so that we don't divide by 0                                                
    sums[sums==0.0]=1.0
    D = sps.spdiags(1.0/sums.flatten(), [0], mtogether.get_shape()[0], mtogether.get_shape()[1], format='csr')
    return D.dot(mtogether)

def random_walk(m_input,t):
    #return m_input.__pow__(t)
    #return np.linalg.matrix_power(m_input,t)
    return m_input.__pow__(t)

def fill_hdf5_with_sparse_by_chunk(mym1,mym2,fname,chunksize):
    start1=0
    end1=0
    n=mym1.shape[0]

    f=h5py.File(fname,'w')
    m1hdf5=f.create_dataset('m1',shape=(n,n),dtype='float')
    m2hdf5=f.create_dataset('m2',shape=(n,n),dtype='float')

    while end1<n:
        end1=np.min([n,(start1+chunksize)])
        print 'start1: '+str(start1)

        if (end1-start1)==1:
            m1hdf5[start1,:]=mym1[start1,:].toarray()
            m2hdf5[start1,:]=mym2[start1,:].toarray()
        else:
            m1hdf5[start1:end1,:]=mym1[start1:end1,:].toarray()
            m2hdf5[start1:end1,:]=mym2[start1:end1,:].toarray()
        start1=end1
    print 'sum of 1'
    print m1hdf5[:,:].sum()
    print m2hdf5[:,:].sum()
    f.close()


def random_walks_by_chunk_get_score_sparse_matrix(mym1,mym2,tmin,tmax,nonzero_total,chunksize):
    scores=[]
    n=mym1.shape[0]
    m1_t=mym1.transpose()
    m2_t=mym2.transpose()

    mat_names[1]='mats'

    for t in range(1,(tmax+1)):
        if t!=1:
            compute_current_matrices(t,mat_names)
        if t>=tmin:
            pass
            #scores.append(1.0*abs_diff_by_chunk_sparse_matrix(t)/nonzero_total)
        print 'done '+str(t)+' '+strftime("%c")
    return scores

def compute_current_matrices(t,mat_names):
    n=mym11.shape[0]

    m12_t=mym12.transpose()
    m22_t=mym22.transpose()

    while end1<n:
        end1=np.min([n,(start1+chunksize)])
        print 'start1: '+str(start1)
      
        if (end1-start1)==1:
            m11_small=mym11[start1,]
            m21_small=mym21[start1,]
        else:
            m11_small=mym11[start1:end1,]
            m21_small=mym21[start1:end1,]
        start1=end1

        start2=0
        end2=0
        while end2<n:
            end2=np.min([n,(start2+chunksize)])
            #print 'start2: '+str(start2)                                                                                                                                   
            #print 'end2: '+str(end2)                                                                                                                                       
            if (end2-start2)==1:
                m12_t_small=m12_t[start2,].transpose()
                m22_t_small=m22_t[start2,].transpose()
            else:
                m12_t_small=m12_t[start2:end2,].transpose()
                m22_t_small=m22_t[start2:end2,].transpose()
            start2=end2

def random_walks_by_chunk_get_score(mym1,mym2,tmin,tmax,nonzero_total,chunksize):
    scores=[]
    hdf5_names={}
    n=mym1.shape[0]
    m1_t=mym1.transpose()
    m2_t=mym2.transpose()

    #write the ms into hdf5s
    #todo: make name more specific
    #print 'filling hdf5 '+strftime("%c")
    hdf5_names[1]='hdf5s'
    fill_hdf5_with_sparse_by_chunk(mym1,mym2,hdf5_names[1],chunksize)

    for t in range(1,(tmax+1)):
        if t!=1:
            hdf5_names[t]='hdf5s_'+str(t)
            #t=1, t=(t-1) and the new t=t that we want to compute
            multiply_by_chunk(hdf5_names[1],hdf5_names[t-1],hdf5_names[t],chunksize)
        if t>=tmin:
            scores.append(1.0*abs_diff_by_chunk(hdf5_names[t],'m1','m2',chunksize)/nonzero_total)
        print 'done '+str(t)+' '+strftime("%c")
    return scores

def get_rss_prop():  # this is quite expensive
    process = psutil.Process(os.getpid())
    return (process.memory_info().rss - process.memory_info().shared) / 10**6

def multiply_by_chunk(hdf5_names_1,hdf5_names_tminus1,hdf5_names_t,chunksize):
    start1=0
    end1=0
    
    f1=h5py.File(hdf5_names_1,'r')
    n=f1['m1'].shape[0]
    f1.close()
    
    ftminus1=h5py.File(hdf5_names_tminus1,'r')
    ft=h5py.File(hdf5_names_t,'w')
    m1t=ft.create_dataset('m1',shape=(n,n),dtype='float')
    m2t=ft.create_dataset('m2',shape=(n,n),dtype='float')
    ftminus1.close()
    ft.close()

    while end1<n:
        end1=np.min([n,(start1+chunksize)])
        print 'start1: '+str(start1)
        print 'memory: '+str(get_rss_prop())

        f1=h5py.File(hdf5_names_1,'r')
        if (end1-start1)==1:
            m11_small=f1['m1'][start1,:]
            m21_small=f1['m2'][start1,:]
        else:
            m11_small=f1['m1'][start1:end1,:]
            m21_small=f1['m2'][start1:end1,:]
        f1.close()

        start2=0
        end2=0
        while end2<n:
            end2=np.min([n,(start2+chunksize)])
            
            ftminus1=h5py.File(hdf5_names_tminus1,'r')
            ft=h5py.File(hdf5_names_t,'r+')
            m1t=ft['m1']
            m2t=ft['m2']
            if (end2-start2)==1:
                m12_small=ftminus1['m1'][:,start2]
                m22_small=ftminus1['m2'][:,start2]
                ftminus1.close()
                if (end1-start1)==1:
                    m1t[start1,start2]+=m11_small[:,:].dot(m12_small[:,:])
                    m2t[start1,start2]+=m21_small[:,:].dot(m22_small[:,:])
                else:
                    m1t[start1:end1,start2]+=m11_small[:,:].dot(m12_small[:,:])
                    m2t[start1:end1,start2]+=m21_small[:,:].dot(m22_small[:,:])
            else:
                m12_small=ftminus1['m1'][:,start2:end2]
                m22_small=ftminus1['m2'][:,start2:end2]
                ftminus1.close()
                if (end1-start1)==1:
                    m1t[start1,start2:end2]+=m11_small[:,:].dot(m12_small[:,:])
                    m2t[start1,start2:end2]+=m21_small[:,:].dot(m22_small[:,:])
                else:
                    m1t[start1:end1,start2:end2]+=m11_small[:,:].dot(m12_small[:,:])
                    m2t[start1:end1,start2:end2]+=m21_small[:,:].dot(m22_small[:,:])
            start2=end2
            ft.close()
            #ftminus1.close()
            
        start1=end1
        #f1.close()
    #f1.close()
    #ftminus1.close()
    #ft.close()

def abs_diff_by_chunk(hdf5_name,m1name,m2name,chunksize):
    start1=0
    end1=0

    f=h5py.File(hdf5_name,'r')
    absdiff=0.0
    n=f[m1name].shape[0]

    print 'm1 sum in diff function'
    print f[m1name][:,:].sum()

    while end1<n:
        end1=np.min([n,(start1+chunksize)])
        print 'start1: '+str(start1)

        if (end1-start1)==1:
            m1_small=f[m1name][start1,:]
            m2_small=f[m2name][start1,:]
        else:
            m1_small=f[m1name][start1:end1,:]
            m2_small=f[m2name][start1:end1,:]
        start1=end1
        
        #print m1_small-m2_small
        print m1_small[:,:].sum()
        absdiff+=abs(m1_small[:,:]-m2_small[:,:]).sum()
        print "===="+str(absdiff)
    return absdiff

def process_by_chunk(mym11,mym12,mym21,mym22,chunksize=1000):
    start1=0
    end1=0
    absdiff=0.0
    n=mym11.shape[0]

    m12_t=mym12.transpose()
    m22_t=mym22.transpose()

    while end1<n:
        end1=np.min([n,(start1+chunksize)])
        print 'start1: '+str(start1)
        #print 'end1: '+str(end1)
        if (end1-start1)==1:
            m11_small=mym11[start1,]
            m21_small=mym21[start1,]
        else:
            m11_small=mym11[start1:end1,]
            m21_small=mym21[start1:end1,]
        start1=end1
    
        start2=0
        end2=0
        while end2<n:
            end2=np.min([n,(start2+chunksize)])
            #print 'start2: '+str(start2)
            #print 'end2: '+str(end2)
            if (end2-start2)==1:
                m12_t_small=m12_t[start2,].transpose()
                m22_t_small=m22_t[start2,].transpose()
            else:
                m12_t_small=m12_t[start2:end2,].transpose()
                m22_t_small=m22_t[start2:end2,].transpose()
            start2=end2
            
            absdiff+=abs(m11_small.dot(m12_t_small)-m21_small.dot(m22_t_small)).sum()
    return absdiff

class DiscoRandomWalks:

    def __init__(self, args):
        self.args = args
    
    def compute_reproducibility(self,m1_csr,m2_csr,args):

        #make symmetric
        m1up=m1_csr
        m1down=m1up.transpose()
        m1down.setdiag(0)
        m1=m1up+m1down

        m2up=m2_csr
        m2down=m2up.transpose()
        m2down.setdiag(0)
        m2=m2up+m2down

        #convert to an actual transition matrix
        if args.transition:
            m1=to_transition(m1)
            m2=to_transition(m2)

        #count nonzero nodes (note that we take the average number of nonzero nodes in the 2 datasets)
	rowsums_1=m1.sum(axis=1)                                                                          
        nonzero_1=[i for i in range(rowsums_1.shape[0]) if rowsums_1[i]>0.0]
	rowsums_2=m2.sum(axis=1)                                                                           
        nonzero_2=[i for i in range(rowsums_2.shape[0]) if rowsums_2[i]>0.0]
	nonzero_total=len(list(set(nonzero_1).union(set(nonzero_2))))
        nonzero_total=0.5*(1.0*len(list(set(nonzero_1)))+1.0*len(list(set(nonzero_2))))

        scores=[]
        big_threshold=30000
        if nonzero_total>big_threshold: #big matrix, process it in chunks
            chunksize=2000
            scores=random_walks_by_chunk_get_score(m1,m2,args.tmin,args.tmax,nonzero_total,chunksize)
        else:
            for t in range(1,args.tmax+1): #range(args.tmin,args.tmax+1):     
                extra_text=' (not included in score calculation)'
                if t==1:
                    rw1=copy.deepcopy(m1)
                    rw2=copy.deepcopy(m2)
                else:
                    rw1=rw1.dot(m1)
                    rw2=rw2.dot(m2)
                if t>=args.tmin:
                    diff=abs(rw1-rw2).sum()#+euclidean(rw1.toarray().flatten(),rw2.toarray().flatten()))
                    scores.append(1.0*float(diff)/float(nonzero_total))
                    extra_text=' | score='+str('{:.3f}'.format(1.0-float(diff)/float(nonzero_total)))
                print 'GenomeDISCO | '+strftime("%c")+' | done t='+str(t)+extra_text
            
                #plot the random walk data
                if not args.concise_analysis:
                    rw1a=rw1.toarray()
                    rw2a=rw2.toarray()
                    combined=np.triu(rw1a)-np.triu(rw2a).T
                    #np.fill_diagonal(combined,0.0)
                    combined=rw1a-rw2a
                    x=mquantiles(abs(combined.flatten()),0.99)
                    plt.matshow(combined,vmin=-x,vmax=x,cmap='bwr')
                    plt.colorbar()
                    #plt.title('Random walk \nt='+str(t))
                    plt.title(args.m1name,fontsize=25,y=1.1,color='red')
                    plt.ylabel(args.m2name,fontsize=25,color='blue')
                    plt.xlabel('random walk iteration '+str(t),fontsize=25)
                    plt.gcf().subplots_adjust(top=0.2)
                    plt.gcf().subplots_adjust(left=0.2)
                    plt.gcf().subplots_adjust(left=0.4)
                    plt.show()
                    fname=args.outdir+'/'+args.outpref+'.'+args.m1name+'.vs.'+args.m2name+'.DiscoRandomWalks.'+str(t)+'.png'
                    plt.savefig(fname)
            
        #compute final score
        ts=range(args.tmin,args.tmax+1)
        denom=len(ts)-1
        if args.tmin==args.tmax:
            auc=scores[0]
        else:
            auc=metrics.auc(range(len(ts)),scores)/denom
        reproducibility=1.0-auc

        #for the report
        reproducibility_text=''
        if not args.concise_analysis:
            reproducibility_text='<br>\n'
            reproducibility_text=reproducibility_text+'<br>\n'
            reproducibility_text=reproducibility_text+'Reproducibility score = '+str(reproducibility)+'\n'
            reproducibility_text=reproducibility_text+'<br>\n'
            rcParams['font.size']= 30
            rcParams['figure.figsize'] = 7,7
            rcParams['xtick.labelsize'] = 20
            rcParams['ytick.labelsize'] = 20
            plt.close("all")
            plt.plot(range(args.tmin,args.tmax+1),scores,'bo-')
            eps=0.02
            for i in range(len(ts)):
                x=range(args.tmin,args.tmax+1)[i]
                y=scores[i]
                plt.text(x,y+eps,str(float("{0:.2f}".format(y))),fontsize=15)
            plt.xlabel('t')
            plt.xticks(range(min(1,args.tmin),args.tmax+2))
            plt.yticks([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
            plt.ylim(0,1.0)
            plt.ylabel('difference/node')
            plt.xlabel('random walk iteration')
            plt.axvline(args.tmax, color='gray', linestyle='dashed', linewidth=2)
            plt.axvline(args.tmin, color='gray', linestyle='dashed', linewidth=2)
            plt.show()
            adj=0.2
            plt.gcf().subplots_adjust(bottom=adj)
            plt.gcf().subplots_adjust(left=adj)
            fname=args.outdir+'/'+args.outpref+'.'+args.m1name+'.vs.'+args.m2name+'.DiscoRandomWalks.Differences.png'
        if not args.concise_analysis:
            plt.savefig(fname)
            reproducibility_text=reproducibility_text+'<img src="'+os.path.basename(fname)+'" width="400" height="400"></td>'+'\n'

        reproducibility_text_rw=''
        if not args.concise_analysis:
            reproducibility_text_rw="<table>"+'\n'
            reproducibility_text_rw=reproducibility_text_rw+'<tr>'+'\n'
            for t in range(args.tmin,args.tmax+1):
                reproducibility_text_rw=reproducibility_text_rw+'<td>'+'\n'
                fname=args.outdir+'/'+args.outpref+'.'+args.m1name+'.vs.'+args.m2name+'.DiscoRandomWalks.'+str(t)+'.png'
                reproducibility_text_rw=reproducibility_text_rw+'<img src="'+os.path.basename(fname)+'" width="400" height="400"></td>'+'\n'
                reproducibility_text_rw=reproducibility_text_rw+'</td>'+'\n'
            reproducibility_text_rw=reproducibility_text_rw+"</table>"+'\n'

        return [reproducibility_text,reproducibility_text_rw],reproducibility,scores

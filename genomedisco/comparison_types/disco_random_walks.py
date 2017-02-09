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
from sklearn import metrics

def random_walk(m_input,t):
    transition=copy.deepcopy(m_input)
    return transition.__pow__(t)

class DiscoRandomWalks:

    def __init__(self, args):
        self.args = args
    
    def compute_reproducibility(self,m1_csr,m2_csr,args):

        #create numpy arrays
        #m1=m1_csr.toarray()
        #m1=m1+m1.T-np.diag(m1.diagonal())
        #m2=m2_csr.toarray()
        #m2=m2+m2.T-np.diag(m2.diagonal())
        #TODO: fix this
        m1=m1_csr+m1_csr.transpose()
        m2=m2_csr+m2_csr.transpose()
	rowsums_1=m1.sum(axis=1)                                                                                                                              
	nonzero_1=[i for i in range(rowsums_1.shape[0]) if rowsums_1[i]>0.0]
	rowsums_2=m2.sum(axis=1)                                                                                                                              
	nonzero_2=[i for i in range(rowsums_2.shape[0]) if rowsums_2[i]>0.0]
	nonzero_total=len(list(set(nonzero_1).union(set(nonzero_2))))

        #perform random walks and make a plot with them
        rw1={}
        rw2={}
        reproducibility_text="<table>"+'\n'
        reproducibility_text=reproducibility_text+'<tr>'+'\n'
        for t in range(args.tmin,args.tmax+1):
            reproducibility_text=reproducibility_text+'<td>'+'\n'
            rw1[t]=random_walk(m1,t).toarray()
            rw2[t]=random_walk(m2,t).toarray()
            #plot these so that we can refer to them in the html
            combined=np.triu(rw1[t])-np.triu(rw2[t]).T
            np.fill_diagonal(combined,0.0)
            x=mquantiles(abs(combined.flatten()),0.99)
            plt.matshow(combined,vmin=-x,vmax=x,cmap='bwr')
            plt.colorbar()
	    plt.title('Random walk \nt='+str(t))
            plt.show()
            fname=args.outdir+'/'+args.outpref+'.'+args.m1name+'.vs.'+args.m2name+'.DiscoRandomWalks.'+str(t)+'.png'
	    print fname
            plt.savefig(fname)
            reproducibility_text=reproducibility_text+'<img src="'+os.path.basename(fname)+'" width="400" height="400"></td>'+'\n'
        reproducibility_text=reproducibility_text+'</tr>'+'\n'
        reproducibility_text=reproducibility_text+"</table>"+'\n'

	scores=[]
	for t in range(args.tmin,args.tmax+1):
            diff=abs(rw1[t]-rw2[t]).sum()
            scores.append(diff)
	ts=range(args.tmin,args.tmax+1)
	denom=len(ts)
	auc=metrics.auc(range(len(ts)),scores)/denom

        return reproducibility_text,auc/nonzero_total

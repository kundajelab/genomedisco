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
    #return m_input.__pow__(t)
    return np.linalg.matrix_power(m_input,t)

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

        nonzero_total=0.5*(1.0*len(list(set(nonzero_1)))+1.0*len(list(set(nonzero_2))))

        m1=m1.toarray()
        m2=m2.toarray()

        #perform random walks and make a plot with them
        scores=[]
        
        reproducibility_text="<table>"+'\n'
        reproducibility_text=reproducibility_text+'<tr>'+'\n'
        for t in range(args.tmin,args.tmax+1):
            print t
            reproducibility_text=reproducibility_text+'<td>'+'\n'
            rw1=random_walk(m1,t)#.toarray()
            rw2=random_walk(m2,t)#.toarray()
            diff=abs(rw1-rw2).sum()
            scores.append(diff/nonzero_total)

            #m1q=mquantiles(rw1[t].flatten(),0.9)
            #m2q=mquantiles(rw2[t].flatten(),0.9)
            #rw1[t][rw1[t]>m1q]=m1q
            #rw2[t][rw2[t]>m2q]=m2q
            #plot these so that we can refer to them in the html
            if not args.concise_analysis:
                combined=np.triu(rw1)-np.triu(rw2).T
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

        if not args.concise_analysis:
            plt.close("all")
            plt.plot(range(args.tmin,args.tmax+1),scores,'bo')
            plt.xlabel('t')
            plt.ylabel('Difference score')
            plt.show()
            fname=args.outdir+'/'+args.outpref+'.'+args.m1name+'.vs.'+args.m2name+'.DiscoRandomWalks.Differences.png'
            plt.savefig(fname)
            reproducibility_text=reproducibility_text+'<img src="'+os.path.basename(fname)+'" width="400" height="200"></td>'+'\n'
        print scores
        ts=range(args.tmin,args.tmax+1)
        denom=len(ts)
	auc=metrics.auc(range(len(ts)),scores)/denom
        reproducibility=1.0-auc
        return reproducibility_text,reproducibility

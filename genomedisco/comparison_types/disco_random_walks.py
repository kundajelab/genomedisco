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
from pylab import rcParams

def random_walk(m_input,t):
    #return m_input.__pow__(t)
    #return np.linalg.matrix_power(m_input,t)
    return m_input.__pow__(t)

class DiscoRandomWalks:

    def __init__(self, args):
        self.args = args
    
    def compute_reproducibility(self,m1_csr,m2_csr,args):
        
        #make symmetric
        m1=m1_csr+m1_csr.transpose()
        m2=m2_csr+m2_csr.transpose()

        #count nonzero nodes (note that we take the average number of nonzero nodes in the 2 datasets)
	rowsums_1=m1.sum(axis=1)                                                                          
        nonzero_1=[i for i in range(rowsums_1.shape[0]) if rowsums_1[i]>0.0]
	rowsums_2=m2.sum(axis=1)                                                                           
        nonzero_2=[i for i in range(rowsums_2.shape[0]) if rowsums_2[i]>0.0]
	nonzero_total=len(list(set(nonzero_1).union(set(nonzero_2))))
        nonzero_total=0.5*(1.0*len(list(set(nonzero_1)))+1.0*len(list(set(nonzero_2))))

        #perform random walks 
        scores=[]        
        for t in range(args.tmin,args.tmax+1):
            rw1=random_walk(m1,t)#.toarray()
            rw2=random_walk(m2,t)#.toarray()
            diff=abs(rw1-rw2).sum()
            scores.append(diff/nonzero_total)

            #plot the random walk data
            if not args.concise_analysis:
                rw1=rw1.toarray()
                rw2=rw2.toarray()
                combined=np.triu(rw1)-np.triu(rw2).T
                np.fill_diagonal(combined,0.0)
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
        auc=metrics.auc(range(len(ts)),scores)/denom
        reproducibility=1.0-auc

        #for the report
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

        return [reproducibility_text,reproducibility_text_rw],reproducibility

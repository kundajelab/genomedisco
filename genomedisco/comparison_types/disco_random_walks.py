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
import psutil
import scipy.sparse as sps
from genomedisco import processing
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MultipleLocator
from decimal import Decimal
#from statsmodels import robust

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

def write_diff_vector_bedfile(diff_vector,nodes,nodes_idx,out_filename):
    out=gzip.open(out_filename,'w')
    for i in range(diff_vector.shape[0]):
        node_name=nodes_idx[i]
        node_dict=nodes[node_name]
        out.write(str(node_dict['chr'])+'\t'+str(node_dict['start'])+'\t'+str(node_dict['end'])+'\t'+node_name+'\t'+str(diff_vector[i][0])+'\n')
    out.close()

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
        if True:
            diff_vector=np.zeros((m1.shape[0],1))
            for t in range(1,args.tmax+1): #range(args.tmin,args.tmax+1):     
                extra_text=' (not included in score calculation)'
                if t==1:
                    rw1=copy.deepcopy(m1)
                    rw2=copy.deepcopy(m2)
                else:
                    rw1=rw1.dot(m1)
                    rw2=rw2.dot(m2)
                if t>=args.tmin:
                    diff_vector+=abs(rw1-rw2).sum(axis=1)
                    diff=abs(rw1-rw2).sum()#+euclidean(rw1.toarray().flatten(),rw2.toarray().flatten()))
                    scores.append(1.0*float(diff)/float(nonzero_total))
                    extra_text=' | score='+str('{:.3f}'.format(1.0-float(diff)/float(nonzero_total)))
                print 'GenomeDISCO | '+strftime("%c")+' | done t='+str(t)+extra_text

        #compute final score
        ts=range(args.tmin,args.tmax+1)
        denom=len(ts)-1
        if args.tmin==args.tmax:
            auc=scores[0]
        else:
            auc=metrics.auc(range(len(ts)),scores)/denom
        reproducibility=1.0-auc

        
        #compute final diff vector
        if not args.concise_analysis:
            #save the difference vector (then the main method will convert it to a bed file)
            final_diff_vector=diff_vector
            if denom>0:
                final_diff_vector=(1.0/denom)*diff_vector
            #write difference vector as a bed file
            diff_vector_file=args.outdir+'/'+args.outpref+'.'+args.m1name+'.vs.'+args.m2name+'.diffScore.bed.gz'
            #nodes,nodes_idx,blacklist_nodes=processing.read_nodes_from_bed(args.node_file,args.blacklist)
            #write_diff_vector_bedfile(final_diff_vector,nodes,nodes_idx,diff_vector_file)
        
        #now, make 1 plot
        if not args.concise_analysis:
            originals=np.triu(m1.toarray())-np.triu(m2.toarray()).T
            rw=np.triu(rw1.toarray())-np.triu(rw2.toarray()).T
            diff_mat=(rw1-rw2).toarray()
            diff_vector=final_diff_vector
            #==================
            figwidth=40
            figheight=8
            #read in resolution
            resolution_file=args.outdir+'/../../../data/metadata/resolution.txt'
            resolution=0.000001*float(open(resolution_file,'r').readlines()[0].split()[0])
            range_originals=[-0.01,0.01]
            range_rw=[-0.01,0.01]
            range_diff_mat=[-0.01,0.01]
            #==================

            #set ticks
            nnodes=originals.shape[0]
            start=0
            ticklist=[]
            ticknames=[]
            ticksize=10.0
            last_tick=-10.0
            current_tick=0.0
            while start<=nnodes:
                if Decimal(str(current_tick))-Decimal(str(last_tick))==Decimal(ticksize) or start==nnodes:
                    ticklist.append(start)
                    ticknames.append(str(1.0*(start*resolution))+' Mb')
                    last_tick+=1.0*ticksize
                current_tick+=1.0*resolution
                start+=1
            
            fig, plots = plt.subplots(1,3)
            fig.set_size_inches(figwidth,figheight)
            
            #original data
            #=============
            colorbar_ticks=[range_originals[0],0,range_originals[1]]
            im1 = plots[0].matshow(originals,vmin=range_originals[0],vmax=range_originals[1],cmap='bwr')
            plots[0].set_title('Original data')
            # Create divider for existing axes instance
            divider = make_axes_locatable(plots[0])
            cax = divider.append_axes("right", size="10%", pad=1.5)
            cbar = plt.colorbar(im1, cax=cax, ticks=MultipleLocator(0.2), format="%.3f",orientation='vertical')
            cax = divider.append_axes("right", size="20%", pad=0.7)
            plt.yticks([])
            plt.xticks([])
            cax.spines['right'].set_visible(False)
            cax.spines['top'].set_visible(False)
            cax.spines['left'].set_visible(False)
            cax.spines['bottom'].set_visible(False)
            cbar.set_ticks(colorbar_ticks)
            cbar.set_ticklabels(colorbar_ticks)
            cbar.ax.tick_params(labelsize=15)
            plots[0].set_xticks([])
            plots[0].set_yticks(ticklist)
            plots[0].set_xticklabels([],size=15)
            plots[0].set_yticklabels(ticknames,size=15)
            plots[0].yaxis.tick_right()
            
            #random walk data
            #================
            colorbar_ticks=[range_rw[0],0,range_rw[1]]
            im1 = plots[1].matshow(rw,vmin=range_rw[0],vmax=range_rw[1],cmap='bwr')
            plots[1].set_title('Smoothed data')
            # Create divider for existing axes instance
            divider = make_axes_locatable(plots[1])
            cax = divider.append_axes("right", size="10%", pad=1.5)
            cbar = plt.colorbar(im1, cax=cax, ticks=MultipleLocator(0.2), format="%.3f",orientation='vertical')
            cax = divider.append_axes("right", size="20%", pad=0.7)
            plt.yticks([])
            plt.xticks([])
            cax.spines['right'].set_visible(False)
            cax.spines['top'].set_visible(False)
            cax.spines['left'].set_visible(False)
            cax.spines['bottom'].set_visible(False)
            cbar.set_ticks(colorbar_ticks)
            cbar.set_ticklabels(colorbar_ticks)
            cbar.ax.tick_params(labelsize=15)
            plots[1].set_xticks([])
            plots[1].set_yticks(ticklist)
            plots[1].set_xticklabels([],size=15)
            plots[1].set_yticklabels(ticknames,size=15)
            plots[1].yaxis.tick_right()
            
            #random walk data differences
            #============================
            colorbar_ticks=[range_diff_mat[0],0,range_diff_mat[1]]
            im1 = plots[2].matshow(diff_mat,vmin=range_diff_mat[0],vmax=range_diff_mat[1],cmap='bwr')
            plots[2].set_title('Diff matrix (smoothed data)')
            # Create divider for existing axes instance
            divider = make_axes_locatable(plots[2])
            cax = divider.append_axes("right", size="20%", pad=1.5)
            plt.plot(diff_vector,range(len(diff_vector)))
            plt.ylim(0,len(diff_vector))
            plt.xlim(0,1.5)
            plt.gca().invert_yaxis()
            plt.yticks(ticklist,[])
            cax = divider.append_axes("right", size="10%", pad=0.7)
            cbar = plt.colorbar(im1, cax=cax, ticks=MultipleLocator(0.2), format="%.3f",orientation='vertical')
            cbar.set_ticks(colorbar_ticks)
            cbar.set_ticklabels(colorbar_ticks)
            cbar.ax.tick_params(labelsize=15)
            plots[2].set_xticks([])
            plots[2].set_yticks(ticklist)
            plots[2].set_xticklabels([],size=15)
            plots[2].set_yticklabels(ticknames,size=15)
            plots[2].yaxis.tick_right()
            fname=args.outdir+'/'+args.outpref+'.'+args.m1name+'.vs.'+args.m2name+'.GenomeDISCO.png'
            plt.savefig(fname)

            
        #for the report
        reproducibility_text=''
        if not args.concise_analysis:
            reproducibility_text='<br>\n'
            reproducibility_text=reproducibility_text+'<br>\n'
            reproducibility_text=reproducibility_text+'Reproducibility score = '+str(reproducibility)+'\n'
            reproducibility_text=reproducibility_text+'<br>\n'
            '''
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
        '''
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

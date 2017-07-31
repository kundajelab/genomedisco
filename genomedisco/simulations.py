
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import pybedtools
import numpy as np
import argparse
import os
import processing
import data_operations
import gzip
import re
import copy
from random import randint

def main():
    parser = argparse.ArgumentParser(description='Code for simulating data for reproducibility analysis.') 
    parser.add_argument('--outdir')
    parser.add_argument('--resolution',type=int)
    parser.add_argument('--maxdist',type=int)
    parser.add_argument('--tadfile')
    parser.add_argument('--nodefile')
    parser.add_argument('--realdatafile')
    parser.add_argument('--tadmeansize',type=int)
    parser.add_argument('--intertadmeandistance',type=int)
    parser.add_argument('--depth',type=int)
    parser.add_argument('--edgenoise')
    parser.add_argument('--nodenoise')
    parser.add_argument('--eps',type=float)
    parser.add_argument('--boundarynoise')
    parser.add_argument('--numsim',type=int)
    parser.add_argument('--distance_dep_data',default='NA')
    parser.add_argument('--distance_dep_tads',default='NA')
    args = parser.parse_args()

    simulate(args)

def read_bed_into_interval(bed):
    regions=[]
    for line in gzip.open(bed):
        items=line.strip().split('\t')
        regions.append(pybedtools.Interval(items[0], int(items[1]), int(items[2])))
    return regions

def simulate(args):

    maxdist_in_nodes=int(1.0*args.maxdist/args.resolution)

    tads=read_bed_into_interval(args.tadfile)
    real_data=read_in_data(args.realdatafile)
    original_tad_boundary_var=0
    original_tadmatrix=tadfile_to_tadmatrix(args.tadfile,original_tad_boundary_var,args.resolution,args.nodefile)
    dd=get_2_distance_dependence_curves(real_data,maxdist_in_nodes,original_tadmatrix)
    dds={}
    if args.distance_dep_data!='NA':
        ddfiles=args.distance_dep_data.split(',')
        ddtads=args.distance_dep_tads.split(',')
        dds={}
        for ddfile_idx in range(len(ddfiles)):
            ddfile=ddfiles[ddfile_idx]
            ddtad=ddtads[ddfile_idx]
            ddtad_matrix=tadfile_to_tadmatrix(ddtad,original_tad_boundary_var,args.resolution,args.nodefile)
            ddata=read_in_data(ddfile)
            print 'done'
            dds[ddfile_idx]=get_2_distance_dependence_curves(ddata,maxdist_in_nodes,ddtad_matrix)
    
    for i in range(args.numsim):
        intro=args.outdir+'/res_'+str(args.resolution)+'.Depth_'+str(args.depth)+'.MaxDist_'+str(args.maxdist)+'.simulatedTADs_mean'+str(args.tadmeansize)+'.S_'+str(i)
        simulatedtadfile=intro+'.gz'
        simulate_tadfile(args.tadmeansize,args.intertadmeandistance,args.resolution,args.nodefile,simulatedtadfile,'simulated')
        simulated_tad_matrix=tadfile_to_tadmatrix(simulatedtadfile,original_tad_boundary_var,args.resolution,args.nodefile)
        #this tad matrix will be used for the noise simulations for node and edge noise
        print i
        for edgenoise in args.edgenoise.split(','):
            for nodenoise in args.nodenoise.split(','):
                for boundarynoise in args.boundarynoise.split(','):
                    ablist=['a','b']
                    for ab_idx in range(len(ablist)):
                        ab=ablist[ab_idx]
                        if boundarynoise>0.0:
                            #simulate a new tad matrix, to which we add the boundary noise
                            simulated_tad_matrix=tadfile_to_tadmatrix(simulatedtadfile,int(boundarynoise),args.resolution,args.nodefile)
                        if len(dds.keys())>0:
                            #we have multiple dist dep curves
                            for ddfile_idx in dds.keys():
                                prob_matrix=get_probability_matrix(simulated_tad_matrix,dds[ddfile_idx],maxdist_in_nodes,float(edgenoise),args.eps,float(nodenoise))
                                sampled_matrix=sample_interactions(prob_matrix,args.depth)
                                ftowrite=intro+'.EN_'+str(edgenoise)+'_eps_'+str(args.eps)+'.NN_'+str(nodenoise)+'.BN_'+str(boundarynoise)+'.'+ab+'.dd_'+str(ddfile_idx)+'.gz'
                                write_matrix(sampled_matrix,ftowrite,args)
                        else:
                            prob_matrix=get_probability_matrix(simulated_tad_matrix,dd,maxdist_in_nodes,float(edgenoise),args.eps,float(nodenoise))
                            sampled_matrix=sample_interactions(prob_matrix,args.depth)
                            ftowrite=intro+'.EN_'+str(edgenoise)+'_eps_'+str(args.eps)+'.NN_'+str(nodenoise)+'.BN_'+str(boundarynoise)+'.'+ab+'.gz'
                            write_matrix(sampled_matrix,ftowrite,args)

def write_matrix(sampled_matrix,fname,args):
    f=gzip.open(fname,'w')
    for i in range(sampled_matrix.shape[0]):
        for j in range(i,sampled_matrix.shape[0]):
            v=sampled_matrix[i,j]
            if v>0.0:
                n1=i*args.resolution+int(args.resolution/2)
                n2=j*args.resolution+int(args.resolution/2)
                f.write(str(n1)+'\t'+str(n2)+'\t'+str(v)+'\n')
    f.close()

def read_in_data(mname_full):
    mat=processing.load_sparse_csr(mname_full).toarray()
    mat=mat + mat.T
    return mat

def get_median_size_of_intervals(intervals,resolution):
    vals=[]
    for i in range(len(intervals)):
        our_interval=intervals[i]
        chromo,start,end=our_interval[0],int(1.0*float(our_interval[1])),int(1.0*float(our_interval[2]))
        vals.append(abs(start-end))
    return np.median(np.array(vals))

#TODO: assumes you provided the correct chromosome for the tadfile and the nodefile
#so, tadfile and nodefile need to refer to the exact same chromosome!
def simulate_tadfile(tad_size,tad_distance,resolution,nodefile,outfile,chrname='simulated'):
    #read in the nodes to learn the dimensions of the TAD matrix
    nodes,nodes_idx=processing.read_nodes_from_bed(nodefile)
    n=len(nodes_idx)
    
    tad_size_n=int(1.0*tad_size/resolution)
    tad_distance_n=int(1.0*tad_distance/resolution)
    
    tad_intervals=[]
    out=gzip.open(outfile,'w')
    
    current_n=0
    while current_n<n:
        #sample a distance between tads
        sampled_distance_between_tads=np.random.poisson(tad_distance_n,1)[0]
        plus_distance=current_n+sampled_distance_between_tads
        if plus_distance>=n:
            break
        
        #sample a tad
        sampled_tad_size=np.random.poisson(tad_size_n,1)[0]
        plus_tad=plus_distance+sampled_tad_size
        if plus_tad>=n:
            break
        this_interval=('chr'+chrname,plus_distance,plus_tad)
        tad_intervals.append(this_interval)
        current_n=plus_tad
        out.write('chr'+chrname+'\t'+str(plus_distance*resolution)+'\t'+str(plus_tad*resolution)+'\n')
        
def tadfile_to_tadmatrix(tadfile,var_boundary_diff_init,resolution,nodefile):
    #read in the nodes to learn the dimensions of the TAD matrix
    nodes,nodes_idx=processing.read_nodes_from_bed(nodefile)
    n=len(nodes_idx)
    
    var_boundary_diff=int(var_boundary_diff_init/resolution)
    
    tads=read_bed_into_interval(tadfile)
    tad_matrix=np.zeros((n,n))
    for tad in tads:
        chromo,start,end=tad[0],int(1.0*float(tad[1])/resolution),int(1.0*float(tad[2])/resolution)
        if var_boundary_diff!=0.0:
            start=min(n,max(0,int(np.random.normal(0,1,1)[0]*var_boundary_diff)+start))
            end=min(n,max(0,int(np.random.normal(0,1,1)[0]*var_boundary_diff)+end))
            if end<start:
                end=start
        for i in range (start,end):
            if i>=n:
                continue
            tad_matrix[i,start:min(n,end)]=1.0
    return tad_matrix

def sample_interactions(prob_matrix,depth):
    new_m=np.zeros(prob_matrix.shape)
    for i in range(new_m.shape[0]):
        for j in range(i,new_m.shape[0]):
            p=prob_matrix[i,j]
            reads=np.random.binomial(depth, p, size=1)[0]
            new_m[i,j]=reads
            new_m[j,i]=reads
    return new_m


def get_probability_matrix(tad_matrix,dd_dict,maxdist,prob_noise,eps,prob_node):
    dd=dd_dict['dd']
    sd=dd_dict['sd']
    prob_m=np.zeros(tad_matrix.shape)
    for i in range(prob_m.shape[0]):
        for j in range(i,min(prob_m.shape[0],i+maxdist+1)):
            d=abs(i-j)
            contact_delta=np.random.randn(1)[0]
            pij=dd['interTAD'][d]+contact_delta*sd['interTAD'][d]
            if tad_matrix[i,j]==1.0 and dd['intraTAD'][d]>0.0:
                pij=dd['intraTAD'][d]+contact_delta*sd['intraTAD'][d]
            pij=max(0.00000000000001,min(0.9999999999,pij))

            noise_addition=0.0
            add_noise=np.random.binomial(1, prob_noise, size=1)[0]
            if add_noise>0.0 and prob_noise!=0.0:
                up_or_down=1
                if np.random.binomial(1, prob_noise, size=1)[0]>0.0:
                    up_or_down=-1
                noise_addition=eps*up_or_down*pij
            pij_final=min(1.0,max(0.0,pij+noise_addition))
            
            prob_m[i,j]=pij_final
            prob_m[j,i]=pij_final
    for i in range(prob_m.shape[0]):
        remove_node=float(np.random.binomial(1, prob_node, size=1)[0])
        if remove_node>0.0:
            prob_m[i,:]=0.0
            prob_m[:,i]=0.0
    total_probs=np.triu(copy.deepcopy(prob_m)).sum()
    return prob_m/total_probs


def get_2_distance_dependence_curves(m,maxdist,tad_matrix):
    n=tad_matrix.shape[0]
    tadmeans=[]
    nontadmeans=[]
    tad_sd=[]
    nontad_sd=[]
    total=0.0
    for d in range(maxdist+1):
        tad_values=[]
        nontad_values=[]
        for i in range(n):
            if (i+d)<n:
                v=m[i,i+d]
            else:
                continue
            total+=v
            is_tad=False
            if tad_matrix[i,i+d]==1.0:
                is_tad=True
            if is_tad:
                tad_values.append(v)
            else:
                nontad_values.append(v)
        tadmeans.append(np.nan_to_num(np.nanmean(np.array(tad_values))))
        nontadmeans.append(np.nan_to_num(np.nanmean(np.array(nontad_values))))
        tad_sd.append(np.nan_to_num(np.nanstd(np.array(tad_values))))
        nontad_sd.append(np.nan_to_num(np.nanstd(np.array(nontad_values))))
    #now, divide by total to get probabilities
    tadprobs=[]
    nontadprobs=[]
    tadprobs_sd=[]
    nontadprobs_sd=[]
    dd={}
    dd['intraTAD']={}
    dd['interTAD']={}
    sd={}
    sd['intraTAD']={}
    sd['interTAD']={}
    for d in range(maxdist+1):
        if d==0:
            dd['intraTAD'][d]=0.0
            dd['interTAD'][d]=0.0
            sd['intraTAD'][d]=0.0
            sd['interTAD'][d]=0.0
        else:
            dd['intraTAD'][d]=tadmeans[d]/total
            dd['interTAD'][d]=nontadmeans[d]/total
            sd['intraTAD'][d]=tad_sd[d]/total
            sd['interTAD'][d]=nontad_sd[d]/total
        tadprobs.append(dd['intraTAD'][d])
        nontadprobs.append(dd['interTAD'][d])
        tadprobs_sd.append(sd['intraTAD'][d])
        nontadprobs_sd.append(sd['interTAD'][d])
        
    plt.plot(np.log(tadprobs)/np.log(10),label='Intra-TAD',color='red')
    plt.plot(np.log(nontadprobs)/np.log(10),label='Inter-TAD',color='blue')
    plt.xlabel('Distance (in nodes)')
    plt.ylabel('Log10(probability of contact)')
    plt.axvline(25,linewidth=1, color = 'lightgray') ### this is 1Mb
    plt.legend()
    plt.show()

    dd_and_sd={}
    dd_and_sd['dd']=dd
    dd_and_sd['sd']=sd

    return dd_and_sd
    


main()

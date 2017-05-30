
import argparse
import logging
import copy
import re
import os
from time import gmtime, strftime
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from pylab import rcParams

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--m1name',default='Matrix15')
    parser.add_argument('--m2name',default='Matrix70')
    parser.add_argument('--out',default='/ifs/scratch/oursu/encode_highres/results/res500000')
    parser.add_argument('--chromos',default='chr21,chr22')
    parser.add_argument('--tmin',type=int,default=3)
    parser.add_argument('--tmax',type=int,default=7)
    args = parser.parse_args()

    header_col='FF0000'
    picsize="200"
    topscores=0.85


    html=open(args.out+'/results/'+args.m1name+'.vs.'+args.m2name+'/genomewide.'+args.m1name+'.vs.'+args.m2name+'.genomedisco.report.html','w')
    
    html.write("<html>"+'\n')
    html.write("<head>"+'\n')
    html.write("<font color=\""+header_col+"\"> <strong>GenomeDISCO | Genomewide report </font></strong>"+'\n')
    html.write("<br>"+'\n')
    html.write("Report generated on "+strftime("%c")+". Code: <a href=\"http://github.com/kundajelab/genomedisco\">http://github.com/kundajelab/genomedisco</a>. Contact: Oana Ursu oursu@stanford.edu"+'\n')
    html.write("<br>"+'\n')
    html.write("<strong>"+args.m1name+" vs "+args.m2name+"</strong>"+'\n')
    html.write("</head>"+'\n')
    html.write("<body>"+'\n')
    html.write("<br>"+'\n')
    html.write("<br>"+'\n')

    #genomewide score
    html.write("<font color=\""+header_col+"\"> <strong>Reproducibility analysis</font></strong>"+'\n')
    html.write("<br>"+'\n')
    score_sum=0.0
    score_num=0
    scores=[]
    chromos=[]
    scoredict={}
    for chromo in args.chromos.split(','):
        score_num+=1
        f=args.out+'/results/'+args.m1name+".vs."+args.m2name+'/genomedisco.'+chromo+'.'+args.m1name+".vs."+args.m2name+'.txt'
        if os.path.isfile(f):
            score=float(open(args.out+'/results/'+args.m1name+".vs."+args.m2name+'/genomedisco.'+chromo+'.'+args.m1name+".vs."+args.m2name+'.txt','r').readlines()[0].split('\t')[2])
            score_sum+=score
            scores.append(score)
            chromos.append(chromo)
            scoredict[chromo]=score
            
    plt.close("all")
    widthfactor=3
    rcParams['figure.figsize'] = 7*widthfactor,7
    rcParams['xtick.labelsize'] = 20
    rcParams['ytick.labelsize'] = 20
    plt.scatter(range(len(chromos)),np.array(scores),s=100)
    plt.xticks(range(len(chromos)),chromos,rotation='vertical')
    plt.xlabel('chromosome',fontsize=30)
    plt.ylim(0.4,1.0)
    plt.axhline(topscores, color='r', linestyle='dashed',linewidth=2,label='threshold for high-quality datasets')
    plt.axhline(np.array(scores).mean(), color='b', linewidth=2,label='genomewide reproducibility for these datasets')
    plt.yticks([0.4,0.5,0.6,0.7,0.8,0.9,1.0])
    plt.ylabel('reproducibility',fontsize=30)
    plt.gcf().subplots_adjust(bottom=0.25)
    plt.gcf().subplots_adjust(left=0.25)
    plt.legend(loc=3,fontsize=25)
    plt.show()
    chrscores=args.out+'/results/'+args.m1name+".vs."+args.m2name+'/'+args.m1name+'.vs.'+args.m2name+'.chrScores.png'
    plt.savefig(chrscores)
    plt.close()

    genomewide_score=score_sum/score_num

    html.write("<td> <strong>What is GenomeDISCO reproducibility?</strong></td>"+'\n')
    html.write("<br>"+'\n')
    html.write("GenomeDISCO (DIfferences between Smoothed COntact maps) computes reproducibility by comparing 2 contact maps at increasing levels of smoothing. The smoothing is done using random walks on graphs. For each dataset, we run random walks of increasing length, and ask what is the probability that we reach node j starting at node i, given a random walk through the network of t steps, or iterations. The key idea is that <strong>if 2 nodes are in contact, then there should be many high-confidence paths connecting them through the network</strong>, even if perhaps the direct edge between them was undersampled. Short random walks provide information about the local network structures, such as loop cliques and subdomains, whereas longer random walks shift the focus toward global structures such as compartments. For each random walk iteration we compare the 2 smoothed contact maps, obtaining an L1 difference in smoothed contact maps."+'\n')
    html.write("<br>"+'\n')
    html.write("In the end, we integrate information across all random walks by computing the area under the curve of L1 differences vs random walk iterations (see difference plot below), resulting in a dissimilarity score between the 2 contact maps of interest. We transform this dissimilarity into a reproducibility score using the formula \"Reproducibility = 1-d\". This yields a reproducibility score between -1 and 1, with higher values indicating similarity (in practice, the range of scores is [0.4,1])."+'\n')
    html.write("<br>"+'\n')
    html.write("GenomeDISCO runs on each chromosome separately. The genomewide score reported below is the average across all chromosomes. <strong> Higher scores are better</strong>.")
    html.write("<br>"+'\n')
    html.write("<br>"+'\n')
    html.write("<font color=\""+header_col+"\"><strong> Your scores</strong></font>"+'\n')
    html.write("<br>"+'\n')
    html.write("<img src=\""+os.path.basename(chrscores)+"\" width=\""+str(int(picsize)*widthfactor)+"\" height=\""+picsize+"\">"+'\n')
    html.write("<br>"+'\n')
    html.write("<br>"+'\n')
    html.write("Reproducibility (genomewide) = "+str(float("{0:.3f}".format(float(genomewide_score))))+'\n')

    if genomewide_score>=topscores:
        outcome='Congratulations! These datasets are highly reproducible.'
    else:
        outcome='These datasets are less reproducible than our empirically defined threshold. This could be due to low sequencing depth, differences in distance dependence curves, noise. Please be cautious with these datasets.'
    html.write("<br>"+'\n')
    html.write("<br>"+'\n')
    html.write("<font color=\"0033FF\"><strong>"+outcome+"</strong></font>"+'\n')
    html.write("<br>"+'\n')
    html.write("<br>"+'\n')
    html.write("<font color=\""+header_col+"\"> <strong>Analysis by chromosome</font></strong>"+'\n')
    html.write("<br>"+'\n')
    html.write("<td> <strong>The difference plot</strong></td>"+'\n')
    html.write(" shows the L1 difference (normalized to the number of nodes) as a function of random walk iteration."+'\n')
    html.write("<br>"+'\n')
    html.write("<td> <strong>The distance dependence plot</strong></td>"+'\n')
    html.write(" shows the probability of contact as a function of linear genomic distance. If 2 datasets have very difference distance dependence curves, their reproducibility will be lower than if they have similar curves."+'\n')
    html.write("<br>"+'\n')
    html.write("<td> <strong>The remaining columns</strong></td>"+'\n')
    html.write(" below show the contact maps after smoothing with random walks. The upper triangular part plotted in red is "+args.m1name+", while the blue is "+args.m2name+". The colorbar shows red values as positive and blue values as negative purely for visualization purposes (in reality all values are positive)."+'\n')
    html.write("<br>"+'\n')
    html.write("<br>"+'\n')
    
    #big table
    html.write("<table border=\"1\" cellpadding=\"10\" cellspacing=\"0\" style=\"border-collapse:collapse;\">"+'\n')
    html.write("<tr>"+'\n')
    html.write("<td> </td>"+'\n')
    html.write("<td> <strong><center>seqdepth</center></strong></td>"+'\n')
    html.write("<td> <strong><center>reproducibility</center></strong></td>"+'\n')
    html.write("<td> <strong><center>difference plot</center></strong></td>"+'\n')
    html.write("<td> <strong><center>distance dependence</center></td>"+'\n')
    for t in range(args.tmin,args.tmax+1):
	html.write("<td><center><strong>Random walk iteration "+str(t)+"</strong></center></td>"+'\n')
    html.write("</tr>"+'\n')

    for chromo in args.chromos.split(','):
	f1=args.out+'/data/edges/'+args.m1name+'/'+args.m1name+'.'+chromo+'.gz'
        f2=args.out+'/data/edges/'+args.m2name+'/'+args.m2name+'.'+chromo+'.gz'
	f=args.out+'/results/'+args.m1name+".vs."+args.m2name+'/genomedisco.'+chromo+'.'+args.m1name+".vs."+args.m2name+'.txt'
        sf=args.out+'/results/'+args.m1name+".vs."+args.m2name+'/genomedisco.'+chromo+'.'+args.m1name+".vs."+args.m2name+'.seqdepth'
        if os.path.isfile(f):
            s1,s2,ssub1,ssub2=open(sf,'r').readlines()[0].strip().split('\t')
            s1=str(float("{0:.2f}".format(float(float(s1)/1000000))))
            s2=str(float("{0:.2f}".format(float(float(s2)/1000000))))
            html.write("<tr>"+'\n')
            html.write("<td> <strong> "+chromo+"</strong></td>"+'\n')
            score=scoredict[chromo]
            html.write("<td> "+args.m1name+": "+str(s1)+' M'+'\n')
            html.write("<br>"+'\n')
            html.write("<br>"+'\n')
            html.write(args.m2name+": "+str(s2)+" M </td>"+'\n')
            html.write("<td> "+str(float("{0:.3f}".format(float(score))))+" </td>"+'\n')
            diffplot=chromo+"."+args.m1name+".vs."+args.m2name+".DiscoRandomWalks.Differences.png"
            html.write("<td> <img src=\""+diffplot+"\" width=\""+picsize+"\" height=\""+picsize+"\"> </td>"+'\n')
            dd=chromo+"."+args.m1name+".vs."+args.m2name+".distDep.png"
            html.write("<td> <img src=\""+dd+"\" width=\""+picsize+"\" height=\""+picsize+"\"> </td>"+'\n')
            for t in range(args.tmin,args.tmax+1):
                pic=chromo+"."+args.m1name+".vs."+args.m2name+".DiscoRandomWalks."+str(t)+".png"
                html.write("<td> <img src=\""+pic+"\" width=\""+picsize+"\" height=\""+picsize+"\"></td>"+'\n')
            html.write("</tr>"+'\n')

    html.write("</table>"+'\n')
    html.write("<br>"+'\n')
    '''
    html.write("<strong>Parameters used</strong>"+'\n')
    html.write("<br>"+'\n')
    for chromo in args.chromos.split(','):
        html.write(chromo+'\n')
        html.write("<br>"+'\n')
        f=args.out+'/results/'+args.m1name+".vs."+args.m2name+'/genomedisco.'+chromo+'.'+args.m1name+".vs."+args.m2name+'.txt'
        if os.path.isfile(f):
            for line in open(args.out+'/results/'+args.m1name+".vs."+args.m2name+'/'+chromo+'.'+args.m1name+".vs."+args.m2name+'.args','r').readlines():
                html.write(line.strip()+' '+'\n')
        html.write("<br>"+'\n')
    '''
    html.write("</body>"+'\n')
    html.write("</html>"+'\n')

main()


import argparse
import subprocess as subp
import os
import gzip
from time import gmtime, strftime
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from pylab import rcParams

def parse_args():
    parser = argparse.ArgumentParser(description='GenomeDISCO main script')

    #individual parsers
    metadata_samples_parser=argparse.ArgumentParser(add_help=False)
    metadata_samples_parser.add_argument('--metadata_samples',required=True,help='required. A file where each row represents a sample, and the entries are "samplename samplefile". Each of these will be processed. Note: each samplename in the file MUST be unique. Each samplefile listed here should follow the format "chr1 n1 chr2 n2 value"')

    metadata_pairs_parser=argparse.ArgumentParser(add_help=False)
    metadata_pairs_parser.add_argument('--metadata_pairs',required=True,help='required')
    
    datatype_parser=argparse.ArgumentParser(add_help=False)
    datatype_parser.add_argument('--datatype',default='hic')

    nodes_parser=argparse.ArgumentParser(add_help=False)
    nodes_parser.add_argument('--nodes',required=True,help='required')

    outdir_parser=argparse.ArgumentParser(add_help=False)
    outdir_parser.add_argument('--outdir',default='NA',required=True,help='required')

    norm_parser=argparse.ArgumentParser(add_help=False)
    norm_parser.add_argument('--norm',default='sqrtvc')

    tmin_parser=argparse.ArgumentParser(add_help=False)
    tmin_parser.add_argument('--tmin',type=int,default=3)
    
    tmax_parser=argparse.ArgumentParser(add_help=False)
    tmax_parser.add_argument('--tmax',type=int,default=7)

    concise_analysis_parser=argparse.ArgumentParser(add_help=False)
    concise_analysis_parser.add_argument('--concise_analysis',action='store_true')

    baits_parser=argparse.ArgumentParser(add_help=False)
    baits_parser.add_argument('--baits',default='NA')

    running_mode_parser=argparse.ArgumentParser(add_help=False)
    running_mode_parser.add_argument('--running_mode',default='NA')

    subparsers = parser.add_subparsers(help='GenomeDISCO help', dest='command')
    subparsers.required = True #http://bugs.python.org/issue9253#msg186387

    #parsers for commands
    all_parser=subparsers.add_parser('run_all',
                            parents=[metadata_samples_parser,metadata_pairs_parser,datatype_parser,nodes_parser,baits_parser,running_mode_parser,outdir_parser,norm_parser,concise_analysis_parser,tmin_parser,tmax_parser],
                            help='Run all steps in the reproducibility analysis with this single command')

    split_parser=subparsers.add_parser('split',
                            parents=[metadata_samples_parser,datatype_parser,nodes_parser,baits_parser,running_mode_parser,outdir_parser],
                            help='(step 1) split files by chromosome')

    reproducibility_parser=subparsers.add_parser('reproducibility',
                            parents=[metadata_pairs_parser,datatype_parser,tmin_parser,tmax_parser,running_mode_parser,outdir_parser,norm_parser,concise_analysis_parser],
                            help='(step 2) compute reproducibility')

    visualize_parser=subparsers.add_parser('visualize',
                            parents=[metadata_pairs_parser,datatype_parser],
                            help='(step 3) create an html report of the results')

    args = vars(parser.parse_args())
    command = args.pop("command", None)
    return command, args

def split_by_chromosome(datatype,metadata_samples,outdir,baits,running_mode,nodes):
    nodes=os.path.abspath(nodes)
    outdir=os.path.abspath(outdir)
    metadata_samples=os.path.abspath(metadata_samples)

    #make the directory structure for the reproducibility analysis
    subp.check_output(['bash','-c','mkdir -p '+outdir+'/scripts'])
    subp.check_output(['bash','-c','mkdir -p '+outdir+'/data/metadata'])
    subp.check_output(['bash','-c','mkdir -p '+outdir+'/data/edges'])
    subp.check_output(['bash','-c','mkdir -p '+outdir+'/data/nodes'])
    subp.check_output(['bash','-c','mkdir -p '+outdir+'/results'])
    
    #make a list of all the chromosomes in the nodes file
    subp.check_output(['bash','-c','zcat -f '+nodes+' | cut -f1 | sort | uniq | awk \'{print "chr"$0}\' | sed \'s/chrchr/chr/g\' | gzip > '+outdir+'/data/metadata/chromosomes.gz'])

    #split the data into chromosomes
    for chromo_line in gzip.open(outdir+'/data/metadata/chromosomes.gz','r').readlines():
        chromo=chromo_line.strip()

        #nodes ===============
        script_nodes_file=outdir+'/scripts/split/nodes/'+chromo+'.nodes.split_files_by_chromosome.sh'
        subp.check_output(['bash','-c','mkdir -p '+os.path.dirname(script_nodes_file)])
        script_nodes=open(script_nodes_file,'w')
        script_nodes.write("#!/bin/sh"+'\n')
        script_nodes.write('source '+bashrc_file+'\n')
        nodefile=outdir+'/data/nodes/nodes.'+chromo+'.gz'

        print 'GenomeDISCO | '+strftime("%c")+' | Splitting nodes '+chromo

        if datatype=='hic':
            script_nodes.write("zcat -f "+nodes+' | awk \'{print "chr"$1"\\t"$2"\\t"$3"\\t"$4"\\tincluded"}\' | sed \'s/chrchr/chr/g\' | awk -v chromosome='+chromo+' \'{if ($1==chromosome) print $0}\' | gzip > '+nodefile+'\n')

        if datatype=='capturec':
            script_nodes.write("zcat -f "+nodes+' | awk \'{print "chr"$1"\\t"$2"\\t"$3"\\t"$4"\\t"$"5}\' | sed \'s/chrchr/chr/g\' | awk -v chromosome='+chromo+' \'{if ($1==chromosome) print $0}\' | gzip > '+nodefile+'.tmp.gz'+'\n')
            script_nodes.write("zcat -f "+baits+' | awk \'{print "chr"$1"\\t"$2"\\t"$3"\\t"$4"\\t"$5}\' | sed \'s/chrchr/chr/g\' | gzip > '+nodefile+'tmp.baits.gz'+'\n')
            script_nodes.write("$mypython "+os.path.abspath("../scripts/annotate_baits.py")+" --nodes "+nodefile+'.tmp.gz  --baits '+nodefile+'tmp.baits.gz --out '+nodefile+'\n')
            script_nodes.write("rm "+nodefile+'.tmp.gz '+nodefile+'tmp.baits.gz'+'\n')

        script_nodes.close()
        run_script(script_nodes_file,running_mode)

        #edges =====================
        for line in open(metadata_samples,'r').readlines():
            items=line.strip().split()
            samplename=items[0]

            print 'GenomeDISCO | '+strftime("%c")+' | Splitting '+samplename+' '+chromo

            samplefile=items[1]
            script_edges_file=outdir+'/scripts/split/'+samplename+'/'+chromo+'.'+samplename+'.split_files_by_chromosome.sh'
            subp.check_output(['bash','-c','mkdir -p '+os.path.dirname(script_edges_file)])
            script_edges=open(script_edges_file,'w')
            script_edges.write("#!/bin/sh"+'\n')
            script_edges.write('source '+bashrc_file+'\n')
            edgefile=outdir+'/data/edges/'+samplename+'/'+samplename+'.'+chromo+'.gz'
            script_edges.write('mkdir -p '+os.path.dirname(edgefile)+'\n')
            script_edges.write('zcat -f '+samplefile+' | awk \'{print "chr"$1"\\t"$2"\\tchr"$3"\\t"$4"\\t"$5}\' | sed \'s/chrchr/chr/g\' | awk -v chromosome='+chromo+' \'{if ($1==chromosome && $3==chromosome) print $2"\\t"$4"\\t"$5}\' | gzip > '+edgefile+'\n')
            script_edges.close()
            run_script(script_edges_file,running_mode)

    print 'GenomeDISCO | '+strftime("%c")+' | ============================='

def run_script(script_name,running_mode):
    subp.check_output(['bash','-c','chmod 755 '+script_name])
    if running_mode=='NA':
        output=subp.check_output(['bash','-c',script_name])
        if output!='':
            print output
    if running_mode=='write_script':
        pass
    if running_mode=='sge':
        memo='10G'
        output=subp.check_output(['bash','-c','qsub -l h_vmem='+memo+' -o '+script_name+'.o -e '+script_name+'.e '+script_name])

def compute_reproducibility(datatype,metadata_pairs,outdir,norm,tmin,tmax,running_mode,concise_analysis):
    outdir=os.path.abspath(outdir)
    metadata_pairs=os.path.abspath(metadata_pairs)

    for chromo_line in gzip.open(outdir+'/data/metadata/chromosomes.gz','r').readlines():
        chromo=chromo_line.strip()
        for line in open(metadata_pairs,'r').readlines():
            items=line.strip().split()
            samplename1,samplename2=items[0],items[1]
            
            print 'GenomeDISCO | '+strftime("%c")+' | Computing reproducibility for '+samplename1+'.vs.'+samplename2+' '+chromo

            script_comparison_file=outdir+'/scripts/reproducibility/'+samplename1+'.vs.'+samplename2+'/'+chromo+'.'+samplename1+'.vs.'+samplename2+'.genomedisco.sh'
            subp.check_output(['bash','-c','mkdir -p '+os.path.dirname(script_comparison_file)])
            script_comparison=open(script_comparison_file,'w')
            script_comparison.write("#!/bin/sh"+'\n')
            script_comparison.write('source '+bashrc_file+'\n')
            f1=outdir+'/data/edges/'+samplename1+'/'+samplename1+'.'+chromo+'.gz'
            f2=outdir+'/data/edges/'+samplename2+'/'+samplename2+'.'+chromo+'.gz'
            nodefile=outdir+'/data/nodes/nodes.'+chromo+'.gz'
            if os.path.isfile(f1) and os.path.getsize(f1)>20:
                if os.path.isfile(f2) and os.path.getsize(f2)>20:
                    concise_analysis_text=''
                    if concise_analysis:
                        concise_analysis_text=' --concise_analysis'
                    outpath=outdir+'/results/'+samplename1+'.vs.'+samplename2
                    subp.check_output(['bash','-c','mkdir -p '+outpath])
                    script_comparison.write("$mypython -W ignore "+os.path.abspath("../genomedisco/genomedisco/compute_reproducibility.py")+" --m1 "+f1+" --m2 "+f2+" --m1name "+samplename1+" --m2name "+samplename2+" --node_file "+nodefile+" --outdir "+outpath+" --outpref "+chromo+" --m_subsample NA --approximation 10000000 --norm "+norm+" --method RandomWalks "+" --tmin "+str(tmin)+" --tmax "+str(tmax)+concise_analysis_text+'\n')
                    script_comparison.close()
                    run_script(script_comparison_file,running_mode)

    print 'GenomeDISCO | '+strftime("%c")+' | ============================='

def visualize(outdir,tmin,tmax,metadata_pairs):
    header_col='FF0000'
    picsize="200"
    topscores=0.85    

    for line in open(metadata_pairs,'r').readlines():
        items=line.strip().split()
        samplename1,samplename2=items[0],items[1]

        print 'GenomeDISCO | '+strftime("%c")+' | Making report for '+samplename1+'.vs.'+samplename2

        html=open(outdir+'/results/'+samplename1+'.vs.'+samplename2+'/report.'+samplename1+'.vs.'+samplename2+'.genomedisco.html','w')
    
        html.write("<html>"+'\n')
        html.write("<head>"+'\n')
        html.write("<font color=\""+header_col+"\"> <strong>GenomeDISCO | Genomewide report </font></strong>"+'\n')
        html.write("<br>"+'\n')
        html.write("Report generated on "+strftime("%c")+'\n')
        html.write("<br>"+'\n')
        html.write("Code: <a href=\"http://github.com/kundajelab/genomedisco\">http://github.com/kundajelab/genomedisco</a>."+'\n')
        html.write("<br>"+'\n')
        html.write("Contact: Oana Ursu oursu@stanford.edu"+'\n')
        html.write("<br>"+'\n')
        html.write("<strong>"+samplename1+" vs "+samplename2+"</strong>"+'\n')
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
        for chromo_line in gzip.open(outdir+'/data/metadata/chromosomes.gz','r').readlines():
            chromo=chromo_line.strip()
            score_num+=1
            f=outdir+'/results/'+samplename1+".vs."+samplename2+'/'+chromo+'.'+samplename1+".vs."+samplename2+'.txt'
            if os.path.isfile(f):
                score=float(open(outdir+'/results/'+samplename1+".vs."+samplename2+'/'+chromo+'.'+samplename1+".vs."+samplename2+'.txt','r').readlines()[0].split('\t')[2])
                score_sum+=score
                scores.append(score)
                chromos.append(chromo)
                scoredict[chromo]=score
            
        plt.close("all")
        widthfactor=3
        rcParams['figure.figsize'] = 10*widthfactor,10
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
        chrscores=outdir+'/results/'+samplename1+".vs."+samplename2+'/'+samplename1+'.vs.'+samplename2+'.chrScores.png'
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
        html.write("<img src=\""+os.path.basename(chrscores)+"\" width=\""+str(int(1.3*int(picsize))*widthfactor)+"\" height=\""+str(1.3*int(picsize))+"\">"+'\n')
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
        html.write(" below show the contact maps after smoothing with random walks. The upper triangular part plotted in red is "+samplename1+", while the blue is "+samplename2+". The colorbar shows red values as positive and blue values as negative purely for visualization purposes (in reality all values are positive)."+'\n')
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
        for t in range(tmin,tmax+1):
            html.write("<td><center><strong>Random walk iteration "+str(t)+"</strong></center></td>"+'\n')
        html.write("</tr>"+'\n')

        for chromo_line in gzip.open(outdir+'/data/metadata/chromosomes.gz','r').readlines():
            chromo=chromo_line.strip()
            f1=outdir+'/data/edges/'+samplename1+'/'+samplename2+'.'+chromo+'.gz'
            f2=outdir+'/data/edges/'+samplename1+'/'+samplename2+'.'+chromo+'.gz'
            f=outdir+'/results/'+samplename1+".vs."+samplename2+'/'+chromo+'.'+samplename1+".vs."+samplename2+'.txt'
            sf=outdir+'/results/'+samplename1+".vs."+samplename2+'/'+chromo+'.'+samplename1+".vs."+samplename2+'.seqdepth'
            if os.path.isfile(f):
                s1,s2,ssub1,ssub2=open(sf,'r').readlines()[0].strip().split('\t')
                s1=str(float("{0:.2f}".format(float(float(s1)/1000000))))
                s2=str(float("{0:.2f}".format(float(float(s2)/1000000))))
                html.write("<tr>"+'\n')
                html.write("<td> <strong> "+chromo+"</strong></td>"+'\n')
                score=scoredict[chromo]
                html.write("<td> "+samplename1+": "+str(s1)+' M'+'\n')
                html.write("<br>"+'\n')
                html.write("<br>"+'\n')
                html.write(samplename2+": "+str(s2)+" M </td>"+'\n')
                html.write("<td> "+str(float("{0:.3f}".format(float(score))))+" </td>"+'\n')
                diffplot=chromo+"."+samplename1+".vs."+samplename2+".DiscoRandomWalks.Differences.png"
                html.write("<td> <img src=\""+diffplot+"\" width=\""+picsize+"\" height=\""+picsize+"\"> </td>"+'\n')
                dd=chromo+"."+samplename1+".vs."+samplename2+".distDep.png"
                html.write("<td> <img src=\""+dd+"\" width=\""+picsize+"\" height=\""+picsize+"\"> </td>"+'\n')
                for t in range(tmin,tmax+1):
                    pic=chromo+"."+samplename1+".vs."+samplename2+".DiscoRandomWalks."+str(t)+".png"
                    html.write("<td> <img src=\""+pic+"\" width=\""+picsize+"\" height=\""+picsize+"\"></td>"+'\n')
                html.write("</tr>"+'\n')

        html.write("</table>"+'\n')
        html.write("<br>"+'\n')
        html.write("</body>"+'\n')
        html.write("</html>"+'\n')

def run_all(datatype,metadata_samples,outdir,baits,running_mode,nodes,metadata_pairs,norm,tmin,tmax,concise_analysis):
    split_by_chromosome(datatype,metadata_samples,outdir,baits,running_mode,nodes)
    compute_reproducibility(datatype,metadata_pairs,outdir,norm,tmin,tmax,running_mode,concise_analysis)
    visualize(outdir,tmin,tmax,metadata_pairs)

def main():
    command_methods = {'split': split_by_chromosome,
                         'reproducibility': compute_reproducibility,
                         'visualize': visualize,
                       'run_all': run_all}
    command, args = parse_args()
    global bashrc_file
    bashrc_file=os.path.abspath("../genomedisco/scripts/bashrc_genomedisco")
    command_methods[command](**args)


if __name__ == "__main__":
    main()

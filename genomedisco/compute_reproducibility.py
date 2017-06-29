
import argparse
import copy
import re
import os
from time import gmtime, strftime

from genomedisco import data_operations, processing, visualization
from genomedisco.comparison_types.disco_random_walks import DiscoRandomWalks
from genomedisco.comparison_types.disco_random_walks_binarized_matrices import DiscoRandomWalks_binarizedMatrices

def main():
    parser = argparse.ArgumentParser(description='Compute reproducibility of 3D genome data')
    parser.add_argument('--datatype',default='hic')
    parser.add_argument('--m1',type=str,default='/srv/gsfs0/projects/kundaje/users/oursu/3d/LA/merged_nodups/processed_data/HIC014.res40000.byChr.chr21.gz')
    parser.add_argument('--m2',type=str,default='/srv/gsfs0/projects/kundaje/users/oursu/3d/LA/merged_nodups/processed_data/HIC001.res40000.byChr.chr21.gz')
    parser.add_argument('--matrix_format',type=str,default='n1n2val',help='c1n1c2n2val')
    parser.add_argument('--node_file',type=str,default='/srv/gsfs0/projects/kundaje/users/oursu/3d/LA/merged_nodups/nodes/Nodes.w40000.chr21.gz')
    parser.add_argument('--remove_diagonal',action='store_true')
    parser.add_argument('--m1name',type=str,default='HIC014')
    parser.add_argument('--m2name',type=str,default='HIC001')
    parser.add_argument('--outdir',type=str,default='OUT')
    parser.add_argument('--outpref',type=str,default='outpref')
    parser.add_argument('--m_subsample',type=str,default='lowest')
    parser.add_argument('--concise_analysis',action='store_true',help='Add this flag to only output the reproducibility score, and not perform the distance dependence analyses.')
    parser.add_argument('--norm',type=str,default='uniform')
    parser.add_argument('--method',type=str,default='RandomWalks')
    parser.add_argument('--tmin',type=int,default=1)
    parser.add_argument('--tmax',type=int,default=3)
    parser.add_argument('--approximation',type=int,default=40000)
    args = parser.parse_args()

    #write_arguments(args)

    os.system('mkdir -p '+args.outdir)

    print "GenomeDISCO | "+strftime("%c")+" | Starting reproducibility analysis"

    nodes,nodes_idx=processing.read_nodes_from_bed(args.node_file)

    print "GenomeDISCO | "+strftime("%c")+" | Loading contact maps"
    m1=processing.construct_csr_matrix_from_data_and_nodes(args.m1,nodes,args.remove_diagonal)
    m2=processing.construct_csr_matrix_from_data_and_nodes(args.m2,nodes,args.remove_diagonal)

    stats={}
    stats[args.m1name]={}
    stats[args.m2name]={}
    stats[args.m1name]['depth']=m1.sum()
    stats[args.m2name]['depth']=m2.sum()

    m1_subsample=copy.deepcopy(m1)
    m2_subsample=copy.deepcopy(m2)
    if args.m_subsample!='NA':
        if args.m_subsample=='lowest':
            if stats[args.m1name]['depth']>=stats[args.m2name]['depth']:
                m_subsample=copy.deepcopy(m2)
            if stats[args.m1name]['depth']<stats[args.m2name]['depth']:
                m_subsample=copy.deepcopy(m1)
        else:
            m_subsample=processing.construct_csr_matrix_from_data_and_nodes(args.m_subsample,nodes,args.remove_diagonal)
        print "GenomeDISCO | "+strftime("%c")+" | Subsampling to the depth of "+args.m_subsample
        print "GenomeDISCO | "+strftime("%c")+" | Subsampling depth = "+str(m_subsample.sum())
        desired_depth=m_subsample.sum()
        #desired_depth=156023
        if m1.sum()>desired_depth:
            m1_subsample=data_operations.subsample_to_depth(m1,desired_depth)
        if m2.sum()>desired_depth:
            m2_subsample=data_operations.subsample_to_depth(m2,desired_depth)

    stats[args.m1name]['subsampled_depth']=m1_subsample.sum()   
    stats[args.m2name]['subsampled_depth']=m2_subsample.sum()

    print "GenomeDISCO | "+strftime("%c")+' | Normalizing with '+args.norm
    m1_norm=data_operations.process_matrix(m1_subsample,args.norm)
    m2_norm=data_operations.process_matrix(m2_subsample,args.norm)

    if not args.concise_analysis:
        #distance dependence analysis
        print "GenomeDISCO | "+strftime("%c")+" | Distance dependence analysis"
        if args.datatype=='hic':
            m1dd=data_operations.get_distance_dep(m1_subsample)
            m2dd=data_operations.get_distance_dep(m2_subsample)
        if args.datatype=='capturec':
            m1dd=data_operations.get_distance_dep_using_nodes_capturec(m1_subsample,nodes,nodes_idx,args.approximation)
            m2dd=data_operations.get_distance_dep_using_nodes_capturec(m2_subsample,nodes,nodes_idx,args.approximation)
        dd_diff=get_dd_diff(m1dd,m2dd)
        visualization.plot_dds([m1dd,m2dd],[args.m1name,args.m2name],args.outdir+'/'+args.outpref+'.'+args.m1name+'.vs.'+args.m2name+'.distDep',args.approximation)

    print "GenomeDISCO | "+strftime("%c")+" | Computing reproducibility score"
    if args.method=='RandomWalks':
        comparer=DiscoRandomWalks(args)
    reproducibility_text,score,scores=comparer.compute_reproducibility(m1_norm,m2_norm,args)

    print "GenomeDISCO | "+strftime("%c")+" | Writing results"
    write_html_report(stats,args,reproducibility_text,score)
    out=open(args.outdir+'/'+args.outpref+'.'+args.m1name+'.vs.'+args.m2name+'.scores.txt','w')
    out.write(args.m1name+'\t'+args.m2name+'\t'+str('{:.3f}'.format(score))+'\n')
    out.close()
    out=open(args.outdir+'/'+args.outpref+'.'+args.m1name+'.vs.'+args.m2name+'.scoresByStep.txt','w')
    t_strings=[]
    score_strings=[]
    t_counter=0
    for t in range(1,(args.tmax+1)):
        if t>=args.tmin:
            score_strings.append(str('{:.3f}'.format(scores[t_counter])))
            t_counter+=1
        else:
            score_strings.append('NA')
        t_strings.append(str(t))
    out.write('#m1'+'\t'+'m2'+'\t'+'\t'.join(t_strings)+'\n')
    out.write(args.m1name+'\t'+args.m2name+'\t'+'\t'.join(score_strings)+'\n')
    out.close()
    out=open(args.outdir+'/'+args.outpref+'.'+args.m1name+'.vs.'+args.m2name+'.datastats.txt','w')
    out.write('#m1name'+'\t'+'m2name'+'\t'+'SeqDepth.m1'+'\t'+'SeqDepth.m2'+'\t'+'SubsampledSeqDepth.m1'+'\t'+'SubsampledSeqDepth.m2'+'\t'+'DistDepDiff'+'\n')
    dd_value='NA'
    if not args.concise_analysis:
        dd_value=str('{:.10f}'.format(dd_diff))
    out.write(args.m1name+'\t'+args.m2name+'\t'+str(stats[args.m1name]['depth'])+'\t'+str(stats[args.m2name]['depth'])+'\t'+str(stats[args.m1name]['subsampled_depth'])+'\t'+str(stats[args.m2name]['subsampled_depth'])+'\t'+dd_value+'\n')
    out.close()
    print "GenomeDISCO | "+'\t'.join(score_strings)
    print "GenomeDISCO | "+strftime("%c")+" | DONE"

def get_dd_diff(m1dd,m2dd):
    d=0.0
    k=set(m1dd.keys()).union(set(m2dd.keys()))
    for key in k:
        if key in m1dd:
            m1val=m1dd[key]
        else:
            m1val=0.0
        if key in m2dd:
            m2val=m2dd[key]
        else:
            m2val=0.0
        d+=abs(m1val-m2val)
    return d

def write_html_report(stats,args,reproducibility_text,score):
    header_col='"009900"'
    #header_col='"#000000"'

    os.system('mkdir -p '+args.outdir)
    out=open(args.outdir+'/'+args.outpref+'.'+args.m1name+'.vs.'+args.m2name+'.report.html','w')

    out.write('<html>'+'\n')
    out.write('<head>'+'\n')
    out.write('<font color='+header_col+'> <strong>Reproducibility report</font></strong>'+'\n')
    out.write('<br>'+'\n')
    out.write('<br>'+'\n')
    #out.write('<strong>'+args.outpref+'. '+'<font color="#FF0000">'+args.m1name+' (red)'+'<font color="#000000">'+' vs '+'<font color="#0000FF">'+args.m2name+' (blue)'+'<font color="#000000">'+'</strong>'+'\n')
    out.write('<strong>'+args.outpref+'. '+args.m1name+' vs '+args.m2name+'</strong>'+'\n')
    out.write('<br>'+'\n')
    out.write('Report generated '+strftime("%Y-%m-%d %H:%M:%S", gmtime())+' by GenomeDISCO (DIfferences in Smoothed COntact maps)'+'\n')
    out.write('<br>'+'\n')
    out.write('</head>'+'\n')

    out.write('<body>'+'\n')
    
    out.write('<br>'+'\n')
    out.write('<font color='+header_col+'> <strong>Sequencing stats</font></strong>'+'\n')
    out.write('<br>'+'\n')
    out.write('<br>'+'\n')
    out.write('<table border="1" cellpadding="10" cellspacing="0" style="border-collapse:collapse;">'+'\n')

    out.write('<tr>')
    out.write('<td> </td>'+'\n')
    out.write('<td> <strong>'+args.m1name+'</strong></td>'+'\n')
    out.write('<td> <strong>'+args.m2name+'</strong></td>'+'\n')
    out.write('</tr>')

    out.write('<tr>')
    out.write('<td> <strong>Sequencing depth</strong></td>'+'\n')
    out.write('<td> '+str(1.0*stats[args.m1name]['depth']/1000000)+' M</td>'+'\n')
    out.write('<td> '+str(1.0*stats[args.m2name]['depth']/1000000)+' M</td>'+'\n')
    out.write('</tr>')

    if args.m_subsample!='NA':
        out.write('<tr>')
        out.write('<td> <strong>Sequencing depth (subsampled) </strong></td>'+'\n')
        out.write('<td> '+str(1.0*stats[args.m1name]['subsampled_depth']/1000000)+' M</td>'+'\n')
        out.write('<td> '+str(1.0*stats[args.m2name]['subsampled_depth']/1000000)+' M</td>'+'\n')
        out.write('</tr>')

    out.write('</table>'+'\n')

    out.write('<br>'+'\n')
    out.write('<font color='+header_col+'><strong>Reproducibility analysis</strong></font>'+'\n')
    out.write(reproducibility_text[0]+'\n')

    if not args.concise_analysis:
        out.write('<br>'+'\n')
        out.write('<font color='+header_col+'><strong>Distance dependence</strong></font>'+'\n')
        out.write('<br>'+'\n')
        out.write('<br>'+'\n')
        out.write('<img src="'+args.outpref+'.'+args.m1name+'.vs.'+args.m2name+'.distDep.png" width="400" height="400"> '+'\n')
        out.write('<br>'+'\n')
        out.write('<br>'+'\n')
        out.write('<font color='+header_col+'><strong>Random walk matrices</strong></font>'+'\n')
        out.write(reproducibility_text[1]+'\n')

    out.write('</body>'+'\n')
    out.write('</html>'+'\n')

def write_arguments(args):
    argout=open(args.outdir+'/'+args.outpref+'.'+args.m1name+'.vs.'+args.m2name+'.args','w')
    argdict=vars(args)
    for k in sorted(argdict):
        argout.write(k+"="+str(argdict[k])+'\n')
    argout.close()

if __name__=="__main__":
    main()


import argparse
import logging
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
    parser.add_argument('--matrix_format',type=str,default='n1n2val')
    parser.add_argument('--node_file',type=str,default='/srv/gsfs0/projects/kundaje/users/oursu/3d/LA/merged_nodups/nodes/Nodes.w40000.chr21.gz')
    parser.add_argument('--remove_diagonal',type=bool, default='True')
    parser.add_argument('--m1name',type=str,default='HIC014')
    parser.add_argument('--m2name',type=str,default='HIC001')
    parser.add_argument('--outdir',type=str,default='OUT')
    parser.add_argument('--outpref',type=str,default='outpref')
    parser.add_argument('--m_subsample',type=str,default='/srv/gsfs0/projects/kundaje/users/oursu/3d/LA/merged_nodups/processed_data/HIC071.res40000.byChr.chr21.gz')
    parser.add_argument('--concise_analysis',action='store_true',help='Add this flag to only output the reproducibility score, and not perform the distance dependence analyses.')
    parser.add_argument('--norm',type=str,default='uniform')
    parser.add_argument('--method',type=str,default='RandomWalks')
    parser.add_argument('--tmin',type=int,default=1)
    parser.add_argument('--tmax',type=int,default=3)
    parser.add_argument('--approximation',type=int,default=40000)
    args = parser.parse_args()

    os.system('mkdir -p'+args.outdir)

    logging.basicConfig(format='%(asctime)s GenomeDisco %(message)s',level=logging.DEBUG)
    logging.info('| main: Starting GenomeDisco')

    logging.info('| main: Reading nodes')
    nodes,nodes_idx=processing.read_nodes_from_bed(args.node_file)

    logging.info('| main: Loading matrices')
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
        logging.info('| main: Subsampling to the depth of '+args.m_subsample)
        m_subsample=processing.construct_csr_matrix_from_data_and_nodes(args.m_subsample,nodes,args.remove_diagonal)
        logging.info('| main: Subsampling depth = '+str(m_subsample.sum()))
        desired_depth=m_subsample.sum()
        if m1.sum()>desired_depth:
            m1_subsample=data_operations.subsample_to_depth(m1,desired_depth)
        if m2.sum()>desired_depth:
            m2_subsample=data_operations.subsample_to_depth(m2,desired_depth)

    stats[args.m1name]['subsampled_depth']=m1_subsample.sum()   
    stats[args.m2name]['subsampled_depth']=m2_subsample.sum()

    logging.info('| main: Normalizing with '+args.norm)
    m1_norm=data_operations.process_matrix(m1_subsample,args.norm)
    m2_norm=data_operations.process_matrix(m2_subsample,args.norm)

    if not args.concise_analysis:
        #distance dependence analysis
        logging.info('| main: Distance dependence analysis')
        if args.datatype=='hic':
            m1dd=data_operations.get_distance_dep(m1_subsample)
            m2dd=data_operations.get_distance_dep(m2_subsample)
        if args.datatype=='capturec':
            m1dd=data_operations.get_distance_dep_using_nodes_capturec(m1_subsample,nodes,nodes_idx,args.approximation)
            m2dd=data_operations.get_distance_dep_using_nodes_capturec(m2_subsample,nodes,nodes_idx,args.approximation)
        visualization.plot_dds([m1dd,m2dd],[args.m1name,args.m2name],args.outdir+'/'+args.outpref+'.distDep',args.approximation)

    logging.info('| main: Reproducibility analysis')
    if args.method=='RandomWalks':
        comparer=DiscoRandomWalks(args)
    if args.method=='RandomWalks_binarizedMatrices':
        comparer=DiscoRandomWalks_binarizedMatrices(args)
    reproducibility_text,score=comparer.compute_reproducibility(m1_norm,m2_norm,args)

    logging.info('| main: Writing report')
    write_html_report(stats,args,reproducibility_text,score)
    out=open(args.outdir+'/'+args.outpref+args.m1name+'.vs.'+args.m2name+'.'+args.method+'.Score.txt','w')
    out.write(args.m1name+'\t'+args.m2name+'\t'+str(score)+'\t'+str(stats[args.m1name]['depth'])+'\t'+str(stats[args.m2name]['depth'])+'\t'+str(stats[args.m1name]['subsampled_depth'])+'\t'+str(stats[args.m2name]['subsampled_depth'])+'\n')
    logging.info('| main: DONE!')


def write_html_report(stats,args,reproducibility_text,score):
    os.system('mkdir -p '+args.outdir)
    out=open(args.outdir+'/'+args.outpref+'.report.html','w')

    out.write('<html>'+'\n')
    out.write('<head>'+'\n')
    out.write('<strong>'+args.m1name+' (red) vs '+args.m2name+' (blue)</strong>'+'\n')
    out.write('<br>'+'\n')
    out.write('Report generated '+strftime("%Y-%m-%d %H:%M:%S", gmtime())+' by GenomeDISCO (DIfferences in Smoothed COntact maps)'+'\n')
    out.write('</head>'+'\n')

    out.write('<body>'+'\n')
    
    out.write('<br>'+'\n')
    out.write('<br>'+'\n')
    out.write('<font color="#a569bd"> <strong>General stats</font></strong>'+'\n')

    out.write('<table border="1" cellpadding="10" cellspacing="0" style="border-collapse:collapse;">'+'\n')

    out.write('<tr>')
    out.write('<td> </td>'+'\n')
    out.write('<td> <strong>'+args.m1name+'</strong></td>'+'\n')
    out.write('<td> <strong>'+args.m2name+'</strong></td>'+'\n')
    out.write('</tr>')

    out.write('<tr>')
    out.write('<td> <strong>Sequencing depth</strong></td>'+'\n')
    out.write('<td> '+str(stats[args.m1name]['depth'])+'</td>'+'\n')
    out.write('<td> '+str(stats[args.m2name]['depth'])+'</td>'+'\n')
    out.write('</tr>')

    out.write('<tr>')
    out.write('<td> <strong>Sequencing depth (subsampled) </strong></td>'+'\n')
    out.write('<td> '+str(stats[args.m1name]['subsampled_depth'])+'</td>'+'\n')
    out.write('<td> '+str(stats[args.m2name]['subsampled_depth'])+'</td>'+'\n')
    out.write('</tr>')

    out.write('</table>'+'\n')

    if not args.concise_analysis:
        out.write('<br>'+'\n')
        out.write('<font color="#a569bd"><strong>Distance dependence</strong></font>'+'\n')
        out.write('<br>'+'\n')
        out.write('<img src="'+args.outpref+'.distDep.png" width="200" height="200"> '+'\n')
        out.write('<br>'+'\n')
    
    out.write('<br>'+'\n')
    out.write('<font color="#a569bd"><strong>Reproducibility analysis</strong></font>'+'\n')

    out.write(reproducibility_text+'\n')
    
    out.write('<br>'+'\n')
    out.write('Reproducibility = '+str(score)+'\n')

    out.write('</body>'+'\n')
    out.write('</html>'+'\n')

main()

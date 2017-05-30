
usage()
{
cat <<EOF
usage: `basename $0` options
Splits a set of nodes and edges into intra-chromosomal interaction matrices.
Contact: oursu@stanford.edu
OPTIONS:
   -h     Show this message and exit
   -t Data type. Can be "hic", "capturec". DEFAULT: hic.
   -i Input metadata. A file where each row represents a sample, and the entries are "samplename samplefile". Each of these will be processed. Note: each samplename in the file MUST be unique. Each samplefile listed here should follow the format "chr1 n1 chr2 n2 value"
   -n Nodes. Bed file with the nodes in this study. The fourth column MUST correspond to the values used for n1 and n2 in the input files. 
   -b If -t is set to capturec, then provide here a file containing just the baits. This is used to determine which nodes to use in the calculation of the distance dependence function. If this is left blank, all nodes are used for the calculation.
   -j Whether the files in -m should be processed in parallel by submitting one job per file. Can be "sge" (using sge for submitting jobs) or "not_parallel". DEFAULT: not_parallel
   -o Output path. DEFAULT: OUT
EOF
}

datatype="hic"
i=''
nodes=''
baits=''
out="OUT"
j="not_parallel"
bashrc=$(dirname "$0")/bashrc_genomedisco


while getopts "ht:i:n:b:j:o:" opt
do
    case $opt in
	h)
            usage; exit;;
	t)
	    datatype=$OPTARG;;
	i)
            i=$OPTARG;;
	n)
	    nodes=$OPTARG;;
	b)
	    baits=$OPTARG;;
	j)
	    j=$OPTARG;;
	o) 
            out=$OPTARG;;
	?)
	    usage
	    exit 1;;
    esac    
done

echo "Your parameters: ======================================================"
echo "Type of data: ${datatype}"
echo "Input metadata file: ${i}"
echo "Bins: ${nodes}"
if [[ ${datatype} == 'capturec' ]];
then
    echo "For capturec, the baits are: ${baits}"
fi
echo "Job submission: ${j}"
echo "Output path: ${out}"
echo "========================================================================"

source ${bashrc}
#======
# Setup
#======
mkdir -p ${out}/scripts
mkdir -p ${out}/data
mkdir -p ${out}/data/edges
mkdir -p ${out}/data/nodes
mkdir -p ${out}/results

#================================
# Split input files into chromosomes
#================================

chromos=$(zcat -f ${nodes} | cut -f1 | sort | uniq | awk '{print "chr"$0}' | sed 's/chrchr/chr/g')

#=========================
#split nodes by chromosome
#=========================

s=${out}/scripts/nodes.split_files_by_chromosome.sh
nodefile=${out}/data/nodes/$(basename ${nodes} | sed 's/[.].gz//g')
echo "#!/bin/sh" > ${s}
echo "source ${bashrc}" >> ${s}
for chromo in ${chromos};
do
    if [[ ${datatype} == 'hic' ]];
    then
	echo "zcat -f ${nodes} | awk '{print \"chr\""'$'"1\"\t\""'$'"2\"\t\""'$'"3\"\t\""'$'"4\"\tincluded\"}' | sed 's/chrchr/chr/g' | awk -v chromosome=${chromo} '{if ("'$'"1==chromosome) print "'$'"0}' | gzip > ${nodefile}.${chromo}.gz" >> ${s}
    elif [[ ${datatype} == 'capturec' ]];
    then
	echo "zcat -f ${nodes} | awk '{print \"chr\""'$'"1\"\t\""'$'"2\"\t\""'$'"3\"\t\""'$'"4\"\t\""'$'"5}' | sed 's/chrchr/chr/g' | awk -v chromosome=${chromo} '{if ("'$'"1==chromosome) print "'$'"0}' | gzip > ${nodefile}.${chromo}.gz_tmp.gz" >> ${s}
	echo "zcat -f ${baits} | awk '{print \"chr\""'$'"1\"\t\""'$'"2\"\t\""'$'"3\"\t\""'$'"4\"\t\""'$'"5}' | sed 's/chrchr/chr/g' | gzip > ${nodefile}.${chromo}.gz_tmp.baits.gz" >> ${s}
	echo "${mypython} ${CODEDIR}/scripts/annotate_baits.py --nodes ${nodefile}.${chromo}.gz_tmp.gz --baits ${nodefile}.${chromo}.gz_tmp.baits.gz --out ${nodefile}.${chromo}.gz" >> ${s}
	echo "rm ${nodefile}.${chromo}.gz_tmp.gz ${nodefile}.${chromo}.gz_tmp.baits.gz" >> ${s}
	echo "echo DONE" >> ${s}
    fi
done

#run
run_code ${s} ${j}

#=========================
#split edges by chromosome
#=========================
while read line
do
    read -a items <<< "$line"
    samplename=${items[0]}
    samplefile=${items[1]}
    s=${out}/scripts/${samplename}.split_files_by_chromosome.sh
    echo "#!/bin/sh" > ${s}
    echo "source ${bashrc}" >> ${s}
    out_perchromo=${out}/data/edges/${samplename}/${samplename}
    if [ ! -d "${out}/data/edges/${samplename}" ]; then
	mkdir -p ${out}/data/edges/${samplename}
    fi
    for chromo in ${chromos};
    do
	echo "zcat -f ${samplefile} | awk '{print \"chr\""'$'"1\"\t\""'$'"2\"\tchr\""'$'"3\"\t\""'$'"4\"\t\""'$'"5}' | sed 's/chrchr/chr/g' | awk -v chromosome=${chromo} '{if ("'$'"1==chromosome && "'$'"3==chromosome) print "'$'"2\"\t\""'$'"4\"\t\""'$'"5}' | gzip > ${out_perchromo}.${chromo}.gz" >> ${s} 
	echo "echo DONE" >> ${s}
    done

    #run
    echo $s
    #run_code ${s} ${j}
    

done < ${i}


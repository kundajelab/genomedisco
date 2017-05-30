
usage()
{
cat <<EOF
usage: `basename $0` options
genomedisco: DIfferences in Smoothed COntact maps
Contact: oursu@stanford.edu
OPTIONS:
   -h     Show this message and exit
   -t Type of data. Can be "hic", "capturec". DEFAULT: hic
   -i Input metadata. A list of samples to be compared. Each line is a comparison, with the format "sample1 sample2". These MUST correspond to the sample names provided in splitByChromosome.sh
   -n Nodes. Bed file with the nodes in this study. The fourth column MUST correspond to the values used for n1 and n2 in the input files.                                 
   -j Whether the files in -m should be processed in parallel by submitting one job per file. Can be "sge" (using sge for submitting jobs) or "not_parallel". DEFAULT: not\
_parallel                                                                                                                                                                  
   -o Output path. DEFAULT: OUT       
   -d Distance bin size (in bp). This is used for plotting the distance dependence curve. DEFAULT: 40000 
   -b Bias normalization method. Can be: "uniform" (no normalization), "sqrtvc", or "coverage_norm". DEFAULT: sqrtvc. 
   -r Reproducibility method to use for this analysis. Can be: "RandomWalks". DEFAULT=RandomWalks
   -s tmin. Minimum steps of random walks. DEFAULT=1
   -e tmax. maximum steps of random walks. DEFAULT=3
EOF
}

datatype="hic"
i=''
nodes=''
j="not_parallel"
out="OUT"
distance_bin="40000"
bashrc=$(dirname "$0")/bashrc_genomedisco
normalization="uniform"
method="RandomWalks"
tmin=1
tmax=3

while getopts "ht:i:n:j:o:d:b:r:s:e:" opt
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
	j)
	    j=$OPTARG;;
	o) 
            out=$OPTARG;;
	d)
	    distance_bin=$OPTARG;;
	b) 
	    normalization=$OPTARG;;
	r)
	    method=$OPTARG;;
	s)
	    tmin=$OPTARG;;
	e)
            tmax=$OPTARG;;
	?)
	    usage
	    exit 1;;
    esac    
done

echo "Your parameters: ======================================================"
echo "Datatype: ${datatype}"
echo "Input metadata: ${i}"
echo "Nodes: ${nodes}"
echo "Output path: ${out}"
echo "Distance bin: ${distance_bin}"
echo "Normalization: ${normalization}"
echo "Method: ${method}"
echo "tmin: ${tmin}"
echo "tmax: ${tmax}"
echo "Job submission: ${j}"
echo "========================================================================"

source ${bashrc}
chromos=$(zcat -f ${nodes} | cut -f1 | sort | uniq | awk '{print "chr"$0}' | sed 's/chrchr/chr/g' | sed 's/chr//g' | sort -n | awk '{print "chr"$1}' | paste -sd,)

#=================================
# Make a html report
#=================================
while read line
do   
    read -a items <<< "$line"
    m1name=${items[0]}
    m2name=${items[1]}
    s=${out}/results/genomewide.${m1name}.vs.${m2name}.report.sh
    echo "source ${bashrc}" > ${s}
    echo "${mypython} ${CODEDIR}/scripts/make_report.py --m1name ${m1name} --m2name ${m2name} --out ${out} --chromos ${chromos} --tmin ${tmin} --tmax ${tmax}" >> ${s}
    run_code ${s} ${j}

done < ${i}


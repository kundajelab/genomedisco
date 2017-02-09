
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
chromos=$(zcat -f ${nodes} | cut -f1 | sort | uniq | awk '{print "chr"$0}' | sed 's/chrchr/chr/g')
nodefile=${out}/data/nodes/$(basename ${nodes} | sed 's/[.].gz//g')

#=================================
# Run genomeDisco
#=================================
while read line
do
    read -a items <<< "$line"
    m1name=${items[0]}
    m2name=${items[1]}
    s=${out}/scripts/${m1name}.vs.${m2name}.run_disco.sh
    outdir=${out}/results/${m1name}.vs.${m2name}
    echo "source ${bashrc}" > ${s}
    for chromo in ${chromos};
    do
	f1=${out}/data/edges/${m1name}/${m1name}.${chromo}.gz
	f2=${out}/data/edges/${m2name}/${m2name}.${chromo}.gz
	mkdir -p ${outdir}
	echo "${mypython} ${CODEDIR}/genomedisco/genomeDisco.py --m1 ${f1} --m2 ${f2} --m1name ${m1name} --m2name ${m2name} --node_file ${nodefile}.${chromo}.gz --outdir ${outdir} --outpref ${chromo} --m_subsample NA --approximation ${distance_bin} --norm ${normalization} --method ${method} --tmin ${tmin} --tmax ${tmax}" >> ${s}
    done
    echo ${s}
done < ${i}

exit
#=================================
# Make an html report
#=================================
m1name=$(basename ${hicfile1} | sed 's/[.]gz//g')
m2name=$(basename ${hicfile2} | sed 's/[.]gz//g')
html=${outdir}/${outpref}/analysis/reproducibility/${m1name}.vs.${m2name}.report.html

echo "<html>" > ${html}
echo "<head>" >> ${html}
echo "<strong>${m1name} (red) vs ${m2name} (blue)</strong>" >> ${html}
echo "</head>" >> ${html}
echo "<body>" >> ${html}
echo "<br>" >> ${html}
echo "<br>" >> ${html}
#genome-wide score
echo "<font color=\"#a569bd\"> <strong>Genome-wide reproducibility</font></strong>" >> ${html}
echo "<br>" >> ${html}
genomewide_score=$(zcat -f ${outdir}/${outpref}/analysis/reproducibility/*${m1name}.vs.${m2name}.${method}.Score.txt | cut -f3 | awk '{ sum += $1 } END { if (NR > 0) print sum / NR }')
genomewide_repro=$(echo ${genomewide_score} | awk '{reproscore=1-$1/2}{print reproscore}') 
echo "Genomewide DISCO score (average abs difference per nonzero node): ${genomewide_score}" >> ${html}
echo "<br>" >> ${html}
echo "Genomewide reproducibility: ${genomewide_repro}" >> ${html}
echo "<br>" >> ${html}
echo "<br>" >> ${html}
#reproducibility score
echo "<font color=\"#a569bd\"> <strong>Reproducibility analysis</font></strong>" >> ${html}
echo "<table border=\"1\" cellpadding=\"10\" cellspacing=\"0\" style=\"border-collapse:collapse;\">" >> ${html}
echo "<tr>" >> ${html}
echo "<td> </td>" >> ${html}
echo "<td> <strong>DISCO distance</strong></td>" >> ${html}
echo "<td> <strong>DISCO reproducibility</strong></td>" >> ${html}
echo "<td> <strong>Distance dependence</td>" >> ${html}
for t in {1..3};
do
    echo "<td> </td>" >> ${html}
done
echo "</tr>" >> ${html}
for chromo in ${chromos};
do
    echo "<tr>" >> ${html}
    echo "<td> <strong> ${chromo}</strong></td>" >> ${html}
    score=$(zcat -f ${outdir}/${outpref}/analysis/reproducibility/${chromo}${m1name}.vs.${m2name}.${method}.Score.txt | cut -f3)
    repro=$(echo ${score} | awk '{reproscore=1-$1/2}{print reproscore}')
    echo "<td> ${score} </td>" >> ${html}
    echo "<td> ${repro} </td>" >> ${html}
    dd=${chromo}.distDep.png
    echo "<td> <img src=\"${dd}\" width=\"400\" height=\"400\"> </td>" >> ${html}
    for t in {1..3};
    do
	pic=${chromo}.${m1name}.vs.${m2name}.DiscoRandomWalks.${t}.png
	echo "<td> <img src=\"${pic}\" width=\"400\" height=\"400\"></td>" >> ${html}
    done
    echo "</tr>" >> ${html}
done
echo "</table>" >> ${html}
echo "<br>" >> ${html}
echo "<br>" >> ${html}
#seq depth
echo "<font color=\"#a569bd\"> <strong>Sequencing depth</font></strong>" >> ${html}
echo "<table border=\"1\" cellpadding=\"10\" cellspacing=\"0\" style=\"border-collapse:collapse;\">" >> ${html}
echo "<tr>" >> ${html}
echo "<td> </td>" >> ${html}
echo "<td> <strong>${m1name}</strong></td>" >> ${html}
echo "<td> <strong>${m2name}</strong></td>" >> ${html}
echo "<td> <strong>${m1name} (subsampled)</strong></td>" >> ${html}
echo "<td> <strong>${m2name} (subsampled)</strong></td>" >> ${html}
echo "</tr>" >> ${html}
for chromo in ${chromos};
do
    echo "<tr>" >> ${html}
    echo "<td> <strong> Sequencing depth (${chromo})</strong></td>" >> ${html}
    m1_depth=$(zcat -f ${outdir}/${outpref}/analysis/reproducibility/${chromo}${m1name}.vs.${m2name}.${method}.Score.txt | cut -f4)
    m2_depth=$(zcat -f ${outdir}/${outpref}/analysis/reproducibility/${chromo}${m1name}.vs.${m2name}.${method}.Score.txt | cut -f5)
    m1_depth_sub=$(zcat -f ${outdir}/${outpref}/analysis/reproducibility/${chromo}${m1name}.vs.${m2name}.${method}.Score.txt | cut -f6)
    m2_depth_sub=$(zcat -f ${outdir}/${outpref}/analysis/reproducibility/${chromo}${m1name}.vs.${m2name}.${method}.Score.txt | cut -f7)
    echo "<td> ${m1_depth} </td>" >> ${html}
    echo "<td> ${m2_depth} </td>" >> ${html}
    echo "<td> ${m1_depth_sub} </td>" >> ${html}
    echo "<td> ${m2_depth_sub} </td>" >> ${html}
    echo "</tr>" >> ${html}
done
echo "</table>" >> ${html}
echo "</body>" >> ${html}
echo "</html>" >> ${html}
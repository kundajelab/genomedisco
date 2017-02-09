
usage()
{
cat <<EOF
usage: `basename $0` options
genomedisco: DIfferences in Smoothed COntact maps
Contact: oursu@stanford.edu
OPTIONS:
   -h     Show this message and exit
   -t Type of data. Can be "hic", "capturec". DEFAULT: hic
   -1 Input file 1. Assumes contacts in the form 'chr1 bin1 chr2 bin2 value'
   -2 Input file 2. Same format as -m2
   -b Bins. Bed file with the bins in this study. The fourth column MUST correspond to the values used for bin1 and bin2 in the input files
   -d Output directory. DEFAULT: "OUT"
   -p Output prefix. DEFAULT: "genomedisco"
   -r Resolution (in bp). This is used for plotting the distance dependence curve, for HiC data. DEFAULT: 40000 
   -n Normalization method. Can be: "uniform" (no normalization), "sqrtvc", or "coverage_norm". DEFAULT: sqrtvc
   -m Method to use for this analysis. Can be: "RandomWalks". DEFAULT=RandomWalks
   -s tmin. Minimum steps of random walks. DEFAULT=1
   -e tmax. maximum steps of random walks. DEFAULT=3
EOF
}

datatype='hic'
m1=''
m2=''
nodes=''
outdir='OUT'
outpref='genomedisco'
resolution='40000'
bashrc=$(dirname "$0")/bashrc_genomedisco
normalization='uniform'
method=''
tmin=1
tmax=3

while getopts "ht:1:2:b:d:p:r:n:m:s:e:" opt
do
    case $opt in
	h)
            usage; exit;;
	t)
	    datatype=$OPTARG;;
	1)
            m1=$OPTARG;;
	2)
            m2=$OPTARG;;
	b)
	    nodes=$OPTARG;;
	d) 
            outdir=$OPTARG;;
	p)
	    outpref=$OPTARG;;
	r) 
	    resolution=$OPTARG;;
	n) 
	    normalization=$OPTARG;;
	m)
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
echo "Input file 1: ${m1}"
echo "Input file 2: ${m2}"
echo "Bins: ${nodes}"
echo "Output directory: ${outdir}"
echo "Output prefix: ${outpref}"
echo "Resolution: ${resolution} bp"
echo "Normalization: ${normalization}"
echo "Method: ${method}"
echo "tmin: ${tmin}"
echo "tmax: ${tmax}"
echo "========================================================================"


source ${bashrc}
#======
# Setup
#======
mkdir -p ${outdir}/${outpref}/scripts
mkdir -p ${outdir}/${outpref}/analysis
mkdir -p ${outdir}/${outpref}/analysis/data
nodefile=${outdir}/${outpref}/analysis/data/nodes

#================================
# Split input files into chromosomes
#================================
s=${outdir}/${outpref}/scripts/split_hicfiles_by_chromosome.sh
out1_perchromo=${outdir}/${outpref}/analysis/data/$(basename ${m1})
out2_perchromo=${outdir}/${outpref}/analysis/data/$(basename ${m2})
chromos=$(zcat -f ${nodes} | cut -f1 | sort | uniq | awk '{print "chr"$0}' | sed 's/chrchr/chr/g')
echo "source ${bashrc}" > ${s}
for chromo in ${chromos};
do
    #echo "zcat -f ${m1} | awk '{print \"chr\""'$'"1\"\t\""'$'"2\"\tchr\""'$'"3\"\t\""'$'"4\"\t\""'$'"5}' | sed 's/chrchr/chr/g' | awk -v chromosome=${chromo} '{if ("'$'"1==chromosome && "'$'"3==chromosome) print "'$'"2\"\t\""'$'"4\"\t\""'$'"5}' | gzip > ${out1_perchromo}.${chromo}.gz" >> ${s}
    #echo "zcat -f ${m2} | awk '{print \"chr\""'$'"1\"\t\""'$'"2\"\tchr\""'$'"3\"\t\""'$'"4\"\t\""'$'"5}' | sed 's/chrchr/chr/g' | awk -v chromosome=${chromo} '{if ("'$'"1==chromosome && "'$'"3==chromosome) print "'$'"2\"\t\""'$'"4\"\t\""'$'"5}' | gzip > ${out2_perchromo}.${chromo}.gz" >> ${s}
    echo "zcat -f ${nodes} | awk '{print \"chr\""'$'"1\"\t\""'$'"2\"\t\""'$'"3\"\t\""'$'"4}' | sed 's/chrchr/chr/g' | awk -v chromosome=${chromo} '{if ("'$'"1==chromosome) print "'$'"0}' | gzip > ${nodefile}.${chromo}.gz" >> ${s}
done
chmod 755 ${s}
#### ${s}

#=================================
# Run genomeDisco
#=================================
s=${outdir}/${outpref}/scripts/run_disco.sh
echo "source ${bashrc}" > ${s}
for chromo in ${chromos};
do
    f1=${outdir}/${outpref}/analysis/data/$(basename ${m1}).${chromo}.gz
    f2=${outdir}/${outpref}/analysis/data/$(basename ${m2}).${chromo}.gz
    m1name=$(basename ${m1} | sed 's/[.]gz//g')
    m2name=$(basename ${m2} | sed 's/[.]gz//g')
    mkdir -p ${outdir}/${outpref}/analysis/reproducibility
    echo "${mypython} ${CODEDIR}/genomedisco/genomeDisco.py --m1 ${f1} --m2 ${f2} --m1name ${m1name} --m2name ${m2name} --node_file ${nodefile}.${chromo}.gz --outdir ${outdir}/${outpref}/analysis/reproducibility --outpref ${chromo} --m_subsample NA --approximation ${resolution} --norm ${normalization} --method ${method} --tmin ${tmin} --tmax ${tmax}" >> ${s}
done
chmod 755 ${s}
echo ${s}

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

hicfile1=$1
hicfile2=$2
chrSizes=$3
outdir=$4
outpref=$5
resolution=$6
bashrc=$7
normalization=$8
method=$9

#============================================================================

source ${bashrc}

#======
# Setup
#======
mkdir -p ${outdir}/${outpref}/scripts
mkdir -p ${outdir}/${outpref}/analysis
mkdir -p ${outdir}/${outpref}/analysis/data

#=========================                                                                                                                                                                           
# Make chromosomal windows                                                                                                                                                                           
#=========================                                                                                                                                                                           

s=${outdir}/${outpref}/scripts/make_windows.sh
nodefile=${outdir}/${outpref}/analysis/data/nodefile.res${resolution}
echo "source ${bashrc}" > ${s}
echo "bedtools makewindows -w ${resolution} -g ${chrSizes} | awk '{print "'$'"1\"\t\""'$'"2\"\t\""'$'"3\"\t\""'$'"2}' | gzip > ${nodefile}.gz" >> ${s}
echo "chromos="'$'"(zcat -f ${nodefile}.gz | cut -f1 | sort | uniq)" >> ${s}
echo "for chromo in "'$'"{chromos};do zcat -f ${nodefile}.gz | grep -w "'$'"{chromo} | gzip > ${nodefile}."'$'"{chromo}.gz;done" >> ${s}
chmod 755 ${s}
${s}
chromos=$(zcat -f ${nodefile}.gz | cut -f1 | sort | uniq | sed 's/\n/ /g')
########remove!!!
chromos="chr18 chr19 chr20 chr21 chr22"

#================================
# Split hic file into chromosomes
#================================
s=${outdir}/${outpref}/scripts/split_hicfiles_by_chromosome.sh
out1_perchromo=${outdir}/${outpref}/analysis/data/$(basename ${hicfile1})
out2_perchromo=${outdir}/${outpref}/analysis/data/$(basename ${hicfile2})
echo "source ${bashrc}" > ${s}
for chromo in ${chromos};
do
    echo "zcat -f ${hicfile1} | awk '{print \"chr\""'$'"1\"\t\""'$'"2\"\tchr\""'$'"3\"\t\""'$'"4\"\t\""'$'"5}' | sed 's/chrchr/chr/g' | awk -v chromosome=${chromo} '{if ("'$'"1==chromosome && "'$'"3==chromosome) print "'$'"2\"\t\""'$'"4\"\t\""'$'"5}' | gzip > ${out1_perchromo}.${chromo}.gz" >> ${s}
    echo "zcat -f ${hicfile2} | awk '{print \"chr\""'$'"1\"\t\""'$'"2\"\tchr\""'$'"3\"\t\""'$'"4\"\t\""'$'"5}' | sed 's/chrchr/chr/g' | awk -v chromosome=${chromo} '{if ("'$'"1==chromosome && "'$'"3==chromosome) print "'$'"2\"\t\""'$'"4\"\t\""'$'"5}' | gzip > ${out2_perchromo}.${chromo}.gz" >> ${s}
done
chmod 755 ${s}
${s}

#=================================
# Run genomeDisco
#=================================
s=${outdir}/${outpref}/scripts/run_disco.sh
echo "source ${bashrc}" > ${s}
for chromo in ${chromos};
do
    f1=${outdir}/${outpref}/analysis/data/$(basename ${hicfile1}).${chromo}.gz
    f2=${outdir}/${outpref}/analysis/data/$(basename ${hicfile2}).${chromo}.gz
    m1name=$(basename ${hicfile1} | sed 's/[.]gz//g')
    m2name=$(basename ${hicfile2} | sed 's/[.]gz//g')
    mkdir -p ${outdir}/${outpref}/analysis/reproducibility
    echo "${mypython} ${CODEDIR}/GenomeDisco/genomeDisco.py --m1 ${f1} --m2 ${f2} --m1name ${m1name} --m2name ${m2name} --node_file ${nodefile}.${chromo}.gz --outdir ${outdir}/${outpref}/analysis/reproducibility --outpref ${chromo} --m_subsample NA --resolution ${resolution} --process_matrix ${normalization} --method ${method} --tmin 1 --tmax 3" >> ${s}
done
chmod 755 ${s}
${s}

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
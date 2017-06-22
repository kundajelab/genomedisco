
CODE=/srv/gsfs0/projects/snyder/oursu/software/git/public_genomedisco/genomedisco
bashrc=${CODE}/scripts/bashrc_genomedisco
step=$1

source ${bashrc}
DATA=/srv/gsfs0/projects/kundaje/users/oursu/3d/encode_highres/AllChrAnon
OUT=/ifs/scratch/oursu/encode_highres/results
mkdir -p ${OUT}
chrSizes=/srv/gsfs0/projects/snyder/oursu/data/hg19.chrom.sizes
pairs=/srv/gsfs0/projects/kundaje/users/oursu/3d/encode_highres/ReproducibilityMatrixPairs.txt

#make node files
#===============
if [[ ${step} == "nodes" ]];
then
    mkdir -p ${DATA}/processed/nodes
    for resolution in 10000 40000 500000;
    do
	bedtools makewindows -g ${chrSizes} -w ${resolution} | awk -v reso=${resolution} '{mid=($2+$2+reso)/2}{print $0"\t"mid}' | gzip > ${DATA}/processed/nodes/Nodes.w${resolution}.bed.gz
    done
fi

#stratify what we need to stratify
if [[ ${step} == "stratify" ]];
then
    echo ${pairs}.all
    cat ${pairs} > ${pairs}.all
    for dmin in 1000000 6000000 10000000;
    do
	zcat -f ${pairs} | grep Stratify | awk -v low=${dmin} '{print $1"_"low" "$2"_"low" "$3}' >> ${pairs}.all
    done
    for m in $(zcat -f ${pairs} | grep Stratify | cut -d " " -f1,2 | sed 's/ /\n/g'| sort | uniq);
    do
	echo ${m}
	dmax=1000000000
	for dmin in 1000000 6000000 10000000;
	do
	    new_m=${m}_${dmin}
	    s=${DATA}/${new_m}.int.bed.gz.sh
	    echo "zcat -f ${DATA}/${m}.int.bed.gz | awk -v low=${dmin} -v high=${dmax} '{dist=("'$'"4-"'$'"2)}{if (dist<0) dist=-dist}{if (dist>low && dist<high) print "'$'"0}' | gzip > ${DATA}/${new_m}.int.bed.gz" > ${s}
	    chmod 755 ${s}
	    #cat ${s}
	    #qsub -o ${s}.o -e ${s}.e ${s}
	done
    done
fi

if [[ ${step} == "figure_out_resolution" ]];
then
    mkdir -p ${DATA}/processed/metadata
    metadata=${DATA}/processed/metadata/metadata
    rm ${metadata}*
    for f in $(ls ${DATA}/*gz);
    do
	res="-1"
	for resolution in 5000 20000 250000;
	do
	    if [[ $(zcat -f ${f} | head -n10000 | cut -f2 | grep "${resolution}$"| wc -l) > 0 ]];
	    then
		res=$(echo "test" | awk -v thisres=${resolution} '{r=2*thisres}{print r}')
	    fi
	done
	echo $f
	dataset=$(basename $f | sed 's/.int.bed.gz//g')
	mkdir -p ${DATA}/processed/${dataset}
	echo "${dataset}delim${f}delim${res}" | sed 's/delim/\t/g' >> ${metadata}.res${res}
    done
fi

if [[ ${step} == "split" ]];
then
    for res in 500000 40000 10000;
    do
	metadata=${DATA}/processed/metadata/metadata
	metadata_samples=${metadata}.res${res}
	nodes=${DATA}/processed/nodes/Nodes.w${res}.bed.gz
	outdir=${OUT}/res${res}
	${mypython} ${CODE}/genomedisco/__main__.py split --metadata_samples ${metadata_samples} --datatype hic --nodes ${nodes} --running_mode sge --outdir ${outdir} 
    done
fi

if [[ ${step} == "score_better" ]];
then
    metadata=${DATA}/processed/metadata/metadata
    for res in 10000 40000; #500000;
    do
	rm ${metadata}.res${res}.pairs
	for datapair in $(zcat -f ${pairs}.all | sed 's/ /delim/g');
        do
	    #echo ${datapair}
            m1=$(echo ${datapair} | sed 's/delim/\t/g' | cut -f1)
            m2=$(echo ${datapair} | sed 's/delim/\t/g' | cut -f2)
	    r=$(echo ${datapair} | sed 's/delim/\t/g' | cut -f3)
	    for m in $(zcat -f ${metadata}.res${res} | cut -f1 | sort | uniq);
	    do
		if [[ ${m} == ${m2} ]];
		then
		    echo "${m1}delim${m2}" | sed 's/delim/\t/g' >> ${metadata}.res${res}.pairs
		fi
	    done
	done
	outdir=${OUT}/res${res}
	nodes=${DATA}/processed/nodes/Nodes.w${res}.bed.gz
	metadata_pairs=${metadata}.res${res}.pairs
	cat ${metadata_pairs} | wc -l
	#${mypython} ${CODE}/genomedisco/__main__.py reproducibility --metadata_pairs ${metadata_pairs} --datatype hic --tmin 3 --tmax 3 --outdir ${outdir} --norm sqrtvc --running_mode sge --concise_analysis
    done
fi

if [[ ${step} == "score_again" ]];
then
    memo=" -l h_vmem=50G -l hostname=scg3* "
    for res in 10000;
    do
	for pair in $(ls ${OUT}/res${res}/scripts/reproducibility/);
	do
	    for chromo in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X;
	    do
		scorefile=${OUT}/res${res}/results/genomedisco/${pair}/chr${chromo}.${pair}.scores.txt
		s=${OUT}/res${res}/scripts/reproducibility/${pair}/chr${chromo}.${pair}.genomedisco.sh
		if [ ! -f ${scorefile} ];
		then
		    echo ${s}
		    qsub ${memo} -o ${s}.o -e ${s}.e ${s}
		fi
	    done
	done
    done
fi

if [[ ${step} == "report" ]];
then
    metadata=${DATA}/processed/metadata/metadata
    for res in 40000 10000 500000;
    do
	outdir=${OUT}/res${res}
        nodes=${DATA}/processed/nodes/Nodes.w${res}.bed.gz
        metadata_pairs=${metadata}.res${res}.pairs
	${mypython} ${CODE}/genomedisco/__main__.py visualize --metadata_pairs ${metadata_pairs} --tmin 3 --tmax 3 --outdir ${outdir} 
	mkdir -p ${outdir}/results/genomedisco/summary
	cat ${outdir}/results/genomedisco/*/genomewide_scores.* | head -n1 > ${outdir}/results/genomedisco/summary/genomedisco.res${res}.txt
	cat ${outdir}/results/genomedisco/*/genomewide_scores.* | grep -v "#" >> ${outdir}/results/genomedisco/summary/genomedisco.res${res}.txt
    done
fi

if [[ ${step} == "score" ]];
then
    metadata=${DATA}/processed/metadata/metadata
    for res in 500000;
    do
	#for datapair in $(zcat -f ${metadata}.res${res} | cut -f1,2 | sed 's/\n/ /g');
	
	for datapair in $(zcat -f /srv/gsfs0/projects/kundaje/users/oursu/3d/encode_highres/ReproducibilityMatrixPairs.txt | sed 's/ /delim/g');
	do
	    echo ${datapair}
	    m1=$(echo ${datapair} | sed 's/delim/\t/g' | cut -f1)
	    m2=$(echo ${datapair} | sed 's/delim/\t/g' | cut -f2)
	    echo "${m1} vs ${m2}"
	    
	done
    done
fi

if [[ ${step} == "test_small" ]];
then
    d=/ifs/scratch/oursu/encode_highres/results
    m1=${d}/res500000/data/edges/Matrix6/Matrix6.chr21.gz
    m2=${d}/res500000/data/edges/Matrix63/Matrix63.chr21.gz
    nodes=${d}/res500000/data/nodes/Nodes.w500000.bed.gz.chr21.gz
    outdir=/srv/gsfs0/projects/kundaje/users/oursu/test/testdisco
    s=/srv/gsfs0/projects/kundaje/users/oursu/test/testdisco_small.sh
    echo "source ${bashrc}" > ${s}
    echo "${mypython} ${CODEDIR}/genomedisco/__main__.py --datatype hic --m1 ${m1} --m2 ${m2} --matrix_format n1n2val --node_file ${nodes} --remove_diagonal --m1name m1name --m2name m2name --outdir ${outdir} --outpref pref --norm sqrtvc --approximation 100000 --m_subsample NA --tmin 3 --tmax 7" >> ${s}
    qsub -o ${s}.o -e ${s}.e ${s}
fi

if [[ ${step} == "test" ]];
then
    d=${OUT}
    m1=${d}/res10000/data/edges/Matrix32/Matrix32.chr1.gz
    m2=${d}/res10000/data/edges/Matrix78/Matrix78.chr1.gz
    nodes=${d}/res10000/data/nodes/Nodes.w10000.bed.gz.chr1.gz
    outdir=/srv/gsfs0/projects/kundaje/users/oursu/test/testdisco
    s=/srv/gsfs0/projects/kundaje/users/oursu/test/testdisco.sh
    echo "source ${bashrc}" > ${s}
    echo "${mypython} ${CODEDIR}/genomedisco/__main__.py --datatype hic --m1 ${m1} --m2 ${m2} --matrix_format n1n2val --node_file ${nodes} --remove_diagonal --m1name m1name --m2name m2name --outdir ${outdir} --outpref pref --norm sqrtvc --approximation 100000 --m_subsample NA --concise_analysis" >> ${s}
    qsub -l h_vmem=100G -o ${s}.o -e ${s}.e ${s}
fi


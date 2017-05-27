
CODE=/oak/stanford/groups/akundaje/oursu/code/genomedisco
bashrc=${CODE}/scripts/bashrc_genomedisco
step=$1

source ${bashrc}
DATA=/oak/stanford/groups/akundaje/oursu/3d/reproducibility/encode/encode_2017-05-25/AllChrAnon
OUT=/oak/stanford/groups/akundaje/oursu/3d/reproducibility/encode/encode_2017-05-25/results
mkdir -p ${OUT}
chrSizes=/oak/stanford/groups/akundaje/oursu/general/hg19.chrom.sizes

#make node files
#===============
if [[ ${step} == "nodes" ]];
then
    rm -r ${DATA}/processed/nodes
    mkdir -p ${DATA}/processed/nodes
    for resolution in 10000 40000 500000;
    do
	bedtools makewindows -g ${chrSizes} -w ${resolution} | awk '{mid=($2+$3)/2}{print $0"\t"mid}' | gzip > ${DATA}/processed/nodes/Nodes.w${resolution}.bed.gz
    done
fi

#split into files by chromosome
#==============================
if [[ ${step} == "split" ]];
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
	dataset=$(basename $f | sed 's/.int.bed.gz//g')
	mkdir -p ${DATA}/processed/${dataset}
	echo "${dataset}delim${f}delim${res}" | sed 's/delim/\t/g' >> ${metadata}.res${res}
    done
    for res in 10000 40000 500000;
    do
	nodes=${DATA}/processed/nodes/Nodes.w${res}.bed.gz
	cmd="${CODEDIR}/scripts/splitByChromosome.sh -t hic -i ${metadata}.res${res} -n ${nodes} -j slurm -o ${OUT}/res${res}"
	echo ${cmd}
	${cmd}
    done
fi

#now, let's try to run genomedisco on chr1 for 10kb resolution, to make sure it works
if [[ ${step} == "test" ]];
then
    d=/oak/stanford/groups/akundaje/oursu/3d/reproducibility/encode/encode_2017-05-25/results
    m1=${d}/res10000/data/edges/Matrix32/Matrix32.chr1.gz
    m2=${d}/res10000/data/edges/Matrix78/Matrix78.chr1.gz
    nodes=${d}/res10000/data/nodes/Nodes.w10000.bed.gz.chr1.gz
    outdir=/oak/stanford/groups/akundaje/oursu/test
    ${mypython} ${CODEDIR}/genomedisco/__main__.py --datatype hic --m1 ${m1} --m2 ${m2} --matrix_format n1n2val --node_file ${nodes} --remove_diagonal --m1name m1name --m2name m2name --outdir ${outdir} --outpref pref --norm sqrtvc --approximation 100000 --concise_analysis
fi



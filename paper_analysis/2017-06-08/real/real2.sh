
code=/srv/gsfs0/projects/snyder/oursu/software/git/public_genomedisco/genomedisco
bashrc=${code}/scripts/bashrc_genomedisco
source ${bashrc}
data=/ifs/scratch/oursu/3d/paper/2017-06-08/LA
metadata_samples=${data}/metadata.samples
metadata_pairs=${data}/metadata.pairs
nodes=/srv/gsfs0/projects/kundaje/users/oursu/3d/LA/merged_nodups/nodes/Nodes.w40000.gz
outdir=${data}/reproducibility
mkdir -p ${data}/datasets

step=$1



if [[ ${step} == "metadata" ]];
then
    #metadata for the pairs
    rm ${metadata_pairs}.tmp
    for data1 in $(zcat -f ${metadata_samples} | cut -f1);
    do
	for data2 in $(zcat -f ${metadata_samples} | cut -f1);
	do
	    echo "${data1}delim${data2}" | sed 's/delim/\t/g' >> ${metadata_pairs}.tmp
	done
    done

    ${mypython} ${code}/scripts/orderpairs.py --file ${metadata_pairs}.tmp --out ${metadata_pairs}.tmp2
    cat ${metadata_pairs}.tmp2 | sort | uniq | awk '{if ($1!=$2) print $0}' > ${metadata_pairs}
    rm ${metadata_pairs}.tmp*
fi

if [[ ${step} == "split" ]];
then
    ${mypython} ${code}/genomedisco/__main__.py split --metadata_samples ${metadata_samples} --datatype hic --nodes ${nodes} --running_mode sge --outdir ${outdir}
fi

if [[ ${step} == "run" ]];
then
    ${mypython} ${code}/genomedisco/__main__.py reproducibility --metadata_pairs ${metadata_pairs} --datatype hic --tmin 1 --tmax 7 --outdir ${outdir} --norm sqrtvc --running_mode sge
fi

if [[ ${step} == "scores" ]];
then
    outscores=${outdir}/plots/disco.scores.txt
    mkdir -p ${outdir}/plots
    zcat -f ${outdir}/results/*/chr18*.scores.txt | awk '{print $1"_"$2"\t"$3}' | sort -k1b,1 > ${outscores}.scores
    zcat -f ${outdir}/results/*/chr18*.datastats.txt | awk '{print $1"_"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' | sort -k1b,1 > ${outscores}.stats
    join -1 1 -2 1 ${outscores}.scores ${outscores}.stats | sed 's/ /\t/g' | sed 's/_/\t/g'  > ${outscores}
    rm ${outscores}.scores ${outscores}.stats 
    echo ${outscores}
fi

if [[ ${step} == 'others' ]];
then
    codebase=/srv/gsfs0/projects/snyder/oursu/software/git/public_genomedisco/genomedisco
    bashrc=${codebase}/paper_analysis/2017-02-05/all_methods/methods_bashrc
    source ${bashrc}
    chrSizes=${codebase}/paper_analysis/2017-02-05/all_methods/chrSizes
    resolution=40000

    for chromosome in $(zcat -f ${outdir}/data/metadata/chromosomes.gz | sed 's/\n/ /g' | sed 's/chr//g');
    do
	echo ${chromosome}
	if [[ ${chromosome} == '1' ]];
	then
            echo ${chromosome}
            bins=${outdir}/data/nodes/nodes.chr${chromosome}.gz
            samples=${metadata_samples}.chr${chromosomes}
	    rm ${samples}
	    for sample in $(zcat -f ${metadata_samples} | cut -f1);
	    do
		echo "${sample}delim${outdir}/data/edges/${sample}/${sample}.chr${chromosome}.gz" | sed 's/delim/\t/g' >> ${samples}
	    done
            parameters=${code}/example_parameters.txt
            action=compute
            ${code}/compute_reproducibility.sh -o ${outdir} -n ${bins} -s ${samples} -p ${metadata_pairs} -a compute -b ${bashrc} -r ${resolution} -m ${parameters} -j sge -c chr${chromosome}
	fi
    done
    
fi

code=/srv/gsfs0/projects/snyder/oursu/software/git/public_genomedisco/genomedisco
bashrc=${code}/scripts/bashrc_genomedisco
source ${bashrc}
data=/ifs/scratch/oursu/3d/paper/2017-05-30/LA
metadata_samples=${data}/metadata.multires.samples
metadata_pairs=${data}/metadata.multires.pairs
#nodes=/srv/gsfs0/projects/kundaje/users/oursu/3d/LA/merged_nodups/nodes/Nodes.w40000.gz
outdir=${data}/multires
chrSizes=${code}/paper_analysis/2017-02-05/all_methods/chrSizes

step=$1

if [[ ${step} == "nodes" ]];
then
    for res in 500000 10000;
    do
	bedtools makewindows -i winnum -w ${res} -s ${res} -g ${chrSizes} | awk '{print $1"\t"$2"\t"$3"\t"$2}' | gzip > ${data}/datasets/nodes.${res}.gz
    done
fi

if [[ ${step} == "metadatas" ]];
then
    rm ${metadata_pairs}
    echo "HIC025delimHIC026" | sed 's/delim/\t/g' >> ${metadata_pairs}
    echo "HIC014delimHIC049" | sed 's/delim/\t/g' >> ${metadata_pairs}
    echo "HIC069delimHIC079" | sed 's/delim/\t/g' >> ${metadata_pairs}
    for res in 500000 10000;
    do
	rm ${metadata_samples}.res${res}
	for dataset in HIC025 HIC026 HIC014 HIC049 HIC069 HIC079;
        do
	    f=${data}/datasets/res${res}/${dataset}.res${res}.gz
	    echo "${dataset}delim${f}" | sed 's/delim/\t/g' >> ${metadata_samples}.res${res}
	done
    done
fi

if [[ ${step} == "multiple_resolutions" ]];
then
    for res in 500000 10000;
    do
	for dataset in HIC025 HIC026 HIC014 HIC049 HIC069 HIC079;
	do
	    reads=$(ls /srv/gsfs0/projects/kundaje/users/oursu/3d/LA/merged_nodups/data/*_${dataset}_merged_nodups.txt.gz)
	    echo ${reads}
	    new_file=${data}/datasets/res${res}/${dataset}.res${res}.gz
	    mkdir -p $(dirname ${new_file})
	    s=${new_file}.sh
	    echo "source /srv/gsfs0/projects/snyder/oursu/software/git/genome_utils/3Dutils/bashrc_3D" > ${s}
	    echo "LA_reads_to_n1n2value_bins.sh ${reads} ${new_file} 30 intra ${res}" >> ${s}
	    chmod 755 ${s}
	    qsub -l h_vmem=100G -o ${s}.o -e ${s}.e ${s}
	done
    done

fi



if [[ ${step} == "split" ]];
then
    for res in 500000 10000;
    do
	nodes=${data}/datasets/nodes.${res}.gz
	${mypython} ${code}/genomedisco/__main__.py split --metadata_samples ${metadata_samples}.res${res} --datatype hic --nodes ${nodes} --running_mode sge --outdir ${outdir}/res${res}
    done
fi

if [[ ${step} == "run" ]];
then
    for res in 500000 10000;
    do
	${mypython} ${code}/genomedisco/__main__.py reproducibility --metadata_pairs ${metadata_pairs} --datatype hic --tmin 1 --tmax 7 --outdir ${outdir}/res${res} --norm sqrtvc --running_mode sge
    done
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
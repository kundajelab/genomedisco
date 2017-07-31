
desired_action=$1

#=========================

CODE=/srv/gsfs0/projects/snyder/oursu/software/git/public_genomedisco/genomedisco
bashrc=${CODE}/scripts/bashrc_genomedisco
source ${bashrc}

#=========================
#make some nodes for 50kb
if [[ "${desired_action}" == "nodes" ]];
then
    for res in 50000;
    do
	bashrc_3Dutils=/srv/gsfs0/projects/snyder/oursu/software/git/genome_utils/3Dutils/bashrc_3D
	source ${bashrc_3Dutils}
	nodes=/ifs/scratch/oursu/data/chr21_datasets/nodes.${res}.chr21.gz
	bedtools makewindows -i winnum -w ${res} -s ${res} -g ${chrSizes} | awk -v reso=${res} '{print $1"\t"$2"\t"$3"\t"$2}' | grep -w "chr21" | gzip > ${nodes}
	ls -lh ${nodes}
    done
fi

mdir=/ifs/scratch/oursu/data/chr21_datasets
nodefile=/ifs/scratch/oursu/data/chr21_datasets/nodes.50000.chr21.gz
mnames="GM12878_combined,HMEC,HUVEC,IMR90,K562,KBM7,NHEK"
dddata=${mdir}/GM12878_combined.chr21.RAWobserved.gz

SIMULATION_DIR=/ifs/scratch/oursu/3d/paper/2017-06-08/simulations
mkdir -p ${SIMULATION_DIR}
EDGE_DIR=${SIMULATION_DIR}/EdgeNoise
NODE_DIR=${SIMULATION_DIR}/NodeNoise
B_DIR=${SIMULATION_DIR}/BoundaryNoise
DD_DIR=${SIMULATION_DIR}/DistanceDependence
NONREP_DIR=${SIMULATION_DIR}/RepNonrep
mkdir -p ${EDGE_DIR}
mkdir -p ${NODE_DIR}
mkdir -p ${B_DIR}
mkdir -p ${DD_DIR}
mkdir -p ${NONREP_DIR}
resolution=50000


if [[ "${desired_action}" == "EdgeNoise" ]];
then
    #Edge noise
    for depth in 10000 100000 1000000 10000000;
    do
	for mname in $(echo ${mnames} | sed 's/,/ /g');
	do
	    matpath=${mdir}/${mname}.chr21.RAWobserved.gz
	    edgenoise="0.0,0.1,0.25,0.5,0.75,0.9"
	    nodenoise=0.0
	    boundarynoise=0
	    for edgenoise in 0.0 0.1 0.25 0.5 0.75 0.9;
	    do
		cmd="${mypython} ${CODEDIR}/genomedisco/simulations_from_real_data.py --outdir ${EDGE_DIR} --resolution ${resolution} --nodes ${nodefile} --matrices ${matpath} --matrix_names ${mname} --distDepData ${dddata} --depth ${depth} --edgenoise ${edgenoise} --nodenoise ${nodenoise} --boundarynoise ${boundarynoise}"
		s=${EDGE_DIR}/script_depth${depth}.${mname}.edgenoise$(echo ${edgenoise} | sed 's/,/_/g').sh
		echo "source ${bashrc}" > ${s}
		echo $cmd >> ${s}
		chmod 755 ${s}
		echo $s
		qsub -l hostname=scg3* -o ${s}.o -e ${s}.e ${s} 
	    done
	done
    done
fi

if [[ "${desired_action}" == "NodeNoise" ]];
then
    #Node noise
    for depth in 10000 100000 1000000 10000000; 
    do
	for mname in $(echo ${mnames} | sed 's/,/ /g');
        do
            matpath=${mdir}/${mname}.chr21.RAWobserved.gz
            edgenoise="0.0"
            boundarynoise=0
            for nodenoise in 0.0 0.1 0.25 0.5 0.75 0.9;
            do
                cmd="${mypython} ${CODEDIR}/genomedisco/simulations_from_real_data.py --outdir ${NODE_DIR} --resolution ${resolution} --nodes ${nodefile} --matrices ${matpath} --matrix_names ${mname} --distDepData ${dddata} --depth ${depth} --edgenoise ${edgenoise} --nodenoise ${nodenoise} --boundarynoise ${boundarynoise}"
		s=${NODE_DIR}/script_depth${depth}.${mname}.nodenoise${nodenoise}.sh
		echo "source ${bashrc}" > ${s}
		echo $cmd >> ${s}
		chmod 755 ${s}
		echo $s
		qsub -l hostname=scg3* -o ${s}.o -e ${s}.e ${s}                                            
            done
        done
    done
fi

if [[ "${desired_action}" == "BoundaryNoise" ]];
then
    #Boundary noise
    for depth in 10000 100000 1000000 10000000;
    do
	for mname in $(echo ${mnames} | sed 's/,/ /g');
        do
            matpath=${mdir}/${mname}.chr21.RAWobserved.gz
	    edgenoise=0.0
	    nodenoise=0.0
	    for boundarynoise in 0 1 2 4 8 16 32;
	    do
		cmd="${mypython} ${CODEDIR}/genomedisco/simulations_from_real_data.py --outdir ${B_DIR} --resolution ${resolution} --nodes ${nodefile} --matrices ${matpath} --matrix_names ${mname} --distDepData ${dddata} --depth ${depth} --edgenoise ${edgenoise} --nodenoise ${nodenoise} --boundarynoise ${boundarynoise}"
		s=${B_DIR}/script_depth${depth}.${mname}.boundarynoise.${boundarynoise}.sh
		echo "source ${bashrc}" > ${s}
		echo $cmd >> ${s}
		chmod 755 ${s}
		echo $s
		qsub -l hostname=scg3* -o ${s}.o -e ${s}.e ${s}
	    done
	done
    done
fi

if [[ "${desired_action}" == "DistDep" ]];
then
    #diff distance dependence
    d1=${mdir}/HUVEC.chr21.RAWobserved.gz
    d2=${mdir}/HMEC.chr21.RAWobserved.gz
    for depth in 10000 100000 1000000 10000000;
    do
	for mname in $(echo ${mnames} | sed 's/,/ /g');
        do
            matpath=${mdir}/${mname}.chr21.RAWobserved.gz
	    edgenoise=0.0
	    nodenoise=0.0
	    boundarynoise=0
	    cmd="${mypython} ${CODEDIR}/genomedisco/simulations_from_real_data.py --outdir ${DD_DIR} --resolution ${resolution} --nodes ${nodefile} --matrices ${matpath} --matrix_names ${mname} --distDepData ${d1},${d2} --depth ${depth} --edgenoise ${edgenoise} --nodenoise ${nodenoise} --boundarynoise ${boundarynoise}"

	    s=${DD_DIR}/script_depth${depth}.ddHIC${HICNUM}.${mname}.sh
	    echo "source ${bashrc}" > ${s}
	    echo $cmd >> ${s}
	    chmod 755 ${s}
	    echo $s
	    qsub -l hostname=scg3* -o ${s}.o -e ${s}.e ${s}
	done
    done
fi

edgenoise_metadata=${EDGE_DIR}/Metadata.samples
edgenoise_metadata_pairs=${EDGE_DIR}/Metadata.pairs
nodenoise_metadata=${NODE_DIR}/Metadata.samples
nodenoise_metadata_pairs=${NODE_DIR}/Metadata.pairs
boundarynoise_metadata=${B_DIR}/Metadata.samples
boundarynoise_metadata_pairs=${B_DIR}/Metadata.pairs
nonrep_metadata=${NONREP_DIR}/Metadata.samples
nonrep_metadata_pairs=${NONREP_DIR}/Metadata.pairs
dd_metadata=${DD_DIR}/Metadata.samples
dd_metadata_pairs=${DD_DIR}/Metadata.pairs


if [[ "${desired_action}" == "Metadata" ]];
then
    rm  ${edgenoise_metadata}
    rm ${edgenoise_metadata_pairs}
    nodenoise=0.0
    boundarynoise=0
    for depth in 10000 100000 1000000 10000000;
    do
	for edgenoise in 0.0 0.1 0.25 0.5 0.75 0.9;
	do
            for mname in $(echo ${mnames} | sed 's/,/ /g');
            do
		d1=Depth_${depth}.${mname}.EN_0.0.NN_${nodenoise}.BN_${boundarynoise}.a.dd_0
		d2=Depth_${depth}.${mname}.EN_${edgenoise}.NN_${nodenoise}.BN_${boundarynoise}.b.dd_0
		echo "${d1}delim${d2}delimchr21" | sed 's/delim/\t/g' >> ${edgenoise_metadata_pairs}
		
		d1=Depth_${depth}.${mname}.EN_0.0.NN_${nodenoise}.BN_${boundarynoise}.b.dd_0
                d2=Depth_${depth}.${mname}.EN_${edgenoise}.NN_${nodenoise}.BN_${boundarynoise}.a.dd_0
                echo "${d1}delim${d2}delimchr21" | sed 's/delim/\t/g' >> ${edgenoise_metadata_pairs}

		for ab in a b;
		do
                    chromo="simulated"
		    dataname_short=Depth_${depth}.${mname}.EN_${edgenoise}.NN_${nodenoise}.BN_${boundarynoise}.${ab}.dd_0
                    echo "${dataname_short}delim${EDGE_DIR}/${dataname_short}.gzdelimNA" | sed 's/delim/\t/g' >> ${edgenoise_metadata}
		done
            done
	done
    done

    rm  ${nodenoise_metadata}
    rm ${nodenoise_metadata_pairs}
    edgenoise=0.0
    boundarynoise=0
    for depth in 10000 100000 1000000 10000000;
    do
	for nodenoise in 0.0 0.1 0.25 0.5 0.75 0.9;
        do
	    for mname in $(echo ${mnames} | sed 's/,/ /g');
            do
		d1=Depth_${depth}.${mname}.EN_${edgenoise}.NN_0.0.BN_${boundarynoise}.a.dd_0
                d2=Depth_${depth}.${mname}.EN_${edgenoise}.NN_${nodenoise}.BN_${boundarynoise}.b.dd_0
                echo "${d1}delim${d2}delimchr21" | sed 's/delim/\t/g' >> ${nodenoise_metadata_pairs}

                d1=Depth_${depth}.${mname}.EN_${edgenoise}.NN_0.0.BN_${boundarynoise}.b.dd_0
                d2=Depth_${depth}.${mname}.EN_${edgenoise}.NN_${nodenoise}.BN_${boundarynoise}.a.dd_0
                echo "${d1}delim${d2}delimchr21" | sed 's/delim/\t/g' >> ${nodenoise_metadata_pairs}

		for ab in a b;
		do
                    chromo="simulated"
		    dataname_short=Depth_${depth}.${mname}.EN_${edgenoise}.NN_${nodenoise}.BN_${boundarynoise}.${ab}.dd_0
                    echo "${dataname_short}delim${NODE_DIR}/${dataname_short}.gzdelimNA" | sed 's/delim/\t/g' >> ${nodenoise_metadata}
		done
	    done
	done
    done

    rm  ${boundarynoise_metadata}
    rm ${boundarynoise_metadata_pairs}
    edgenoise=0.0
    nodenoise=0.0
    for depth in 10000 100000 1000000 10000000;
    do
	for boundarynoise in 0 1 2 4 8 16 32;
	do
	    for mname in $(echo ${mnames} | sed 's/,/ /g');
            do
		d1=Depth_${depth}.${mname}.EN_${edgenoise}.NN_${nodenoise}.BN_0.a.dd_0
                d2=Depth_${depth}.${mname}.EN_${edgenoise}.NN_${nodenoise}.BN_${boundarynoise}.b.dd_0
                echo "${d1}delim${d2}delimchr21" | sed 's/delim/\t/g' >> ${boundarynoise_metadata_pairs}

                d1=Depth_${depth}.${mname}.EN_${edgenoise}.NN_${nodenoise}.BN_0.b.dd_0
                d2=Depth_${depth}.${mname}.EN_${edgenoise}.NN_${nodenoise}.BN_${boundarynoise}.a.dd_0
                echo "${d1}delim${d2}delimchr21" | sed 's/delim/\t/g' >> ${boundarynoise_metadata_pairs}

		for ab in a b;
		do
		    chromo="simulated"
                    dataname_short=Depth_${depth}.${mname}.EN_${edgenoise}.NN_${nodenoise}.BN_${boundarynoise}.${ab}.dd_0
                    echo "${dataname_short}delim${B_DIR}/${dataname_short}.gzdelimNA" | sed 's/delim/\t/g' >> ${boundarynoise_metadata}
                done
	    done
	done
    done

    rm  ${nonrep_metadata}
    rm ${nonrep_metadata_pairs}
    boundarynoise=0
    nodenoise=0.0
    edgenoise=0.0
    for depth in 10000 100000 1000000 10000000;
    do
	for edgenoise in 0.0 0.1 0.25 0.5 0.75 0.9;
	do
	    for mname1 in $(echo ${mnames} | sed 's/,/ /g');
            do
		for  ab1 in a b;
		do 
		    dataname_short1=Depth_${depth}.${mname1}.EN_0.0.NN_${edgenoise}.BN_${boundarynoise}.${ab1}.dd_0
		    for mname2 in $(echo ${mnames} | sed 's/,/ /g');
		    do
			for ab2 in a b;
			do
			    dataname_short2=Depth_${depth}.${mname2}.EN_${edgenoise}.NN_${nodenoise}.BN_${boundarynoise}.${ab2}.dd_0
			    echo "${dataname_short1}delim${dataname_short2}delimsimulated" | sed 's/delim/\t/g' >> ${nonrep_metadata_pairs}.many
			done
		    done
		done
	    done
	done
    done
    ${mypython} ${CODEDIR}/scripts/orderpairs.py --file ${nonrep_metadata_pairs}.many --out ${nonrep_metadata_pairs}.many2
    cat ${nonrep_metadata_pairs}.many2 | awk '{if ($1!=$2) print $0}' | sort | uniq > ${nonrep_metadata_pairs}
    cp ${edgenoise_metadata} ${nonrep_metadata}
    rm ${nonrep_metadata_pairs}.many*

    rm  ${dd_metadata}
    rm ${dd_metadata_pairs}.many
    rm ${dd_metadata_pairs}
    boundarynoise=0
    nodenoise=0.0
    edgenoise=0.0
    for depth in 10000 100000 1000000 10000000;
    do
	for mname1 in $(echo ${mnames} | sed 's/,/ /g');
	do
	    for  dd1 in 0 1;
	    do
                for  ab1 in a b;
                do
                    dataname_short1=Depth_${depth}.${mname1}.EN_${edgenoise}.eps_0.9.NN_${nodenoise}.BN_${boundarynoise}.${ab1}.dd_${dd1}
                    for mname2 in ${mname1};
                    do
			for dd2 in 0 1;
			do
                            for ab2 in a b;
                            do
				dataname_short2=Depth_${depth}.${mname2}.EN_${edgenoise}.eps_0.9.NN_${nodenoise}.BN_${boundarynoise}.${ab2}.dd_${dd2}
				echo "${dataname_short1}delim${dataname_short2}delimsimulated" | sed 's/delim/\t/g' >> ${dd_metadata_pairs}.many
			    done
			done
		    done
		done
            done
	done
    done
    cat ${dd_metadata_pairs}.many | awk '{if ($1!=$2) print $0}' | sort | uniq > ${dd_metadata_pairs}
    ${mypython} ${CODEDIR}/scripts/orderpairs.py --file ${dd_metadata_pairs}.many --out ${dd_metadata_pairs}.many2
    cat ${dd_metadata_pairs}.many2 | awk '{if ($1!=$2) print $0}' | sort | uniq > ${dd_metadata_pairs}
    rm ${dd_metadata_pairs}.many*

    for depth in 10000 100000 1000000 10000000;
    do
	for  dd1 in 0 1;
	do
	    for mname in $(echo ${mnames} | sed 's/,/ /g');
            do
                for  ab1 in a b;
                do
                    dataname_short1=Depth_${depth}.${mname}.EN_${edgenoise}.eps_0.9.NN_${nodenoise}.BN_${boundarynoise}.${ab1}.dd_${dd1}
		    chromo="simulated"
		    echo "${dataname_short1}delim${DD_DIR}/${dataname_short1}.gzdelimNA" | sed 's/delim/\t/g' >> ${dd_metadata}
		done
	    done
	done
    done
fi

if [[ "${desired_action}" == "automatic" ]];
then
    for spot in NodeNoise BoundaryNoise; #DistanceDependence BoundaryNoise RepNonrep NodeNoise EdgeNoise;
    do
	chromosome=21
	metadata_pairs=${SIMULATION_DIR}/${spot}/Metadata.pairs 
	samples=${SIMULATION_DIR}/${spot}/Metadata.samples
	outdir=${SIMULATION_DIR}/${spot}
	resolution=40000
	bashrc=${CODEDIR}/paper_analysis/2017-06-08/all_methods/methods_bashrc
	nodes=${nodefile}
	mkdir -p ${outdir}

	cat ${metadata_pairs}
	cat ${samples}
	parameters=${CODEDIR}/paper_analysis/2017-06-08/all_methods/example_parameters.txt
        ${CODEDIR}/paper_analysis/2017-06-08/all_methods/compute_reproducibility.sh -o ${outdir} -n ${nodes} -s ${samples} -p ${metadata_pairs} -a compute -b ${bashrc} -r ${resolution} -m ${parameters} -j sge -c chr${chromosome}
    done
fi


if [[ "${desired_action}" == "scores" ]];
then
    for noisedirname in ${EDGE_DIR} ${NODE_DIR} ${B_DIR} ${DD_DIR} ${NONREP_DIR};
    do
	for method in genomedisco hicrep hic-spector;
	do
	    summary=${noisedirname}/${method}.results.txt
	    rm ${summary} 
	    echo "====="
	    echo ${summary}
	    cat ${summary}
	    addon=""
	    if [[ ${method} == 'genomedisco' ]];
	    then
		addon="scores"
	    fi
	    cat  ${noisedirname}/results/${method}/*/*${addon}.txt | cut -f1,2,3 | sort -k3 -n |sed 's/Depth_//g' | sed 's/[.]G/\tG/g' | sed 's/[.]K/\tK/g' | sed 's/[.]I/\tI/g' | sed 's/[.]H/\tH/g' | sed 's/[.]NH/\tNH/g' | sed 's/[.]EN_/\t/g' | sed 's/[.]NN_/\t/g' | sed 's/[.]eps_0[.]9//g' | sed 's/[.]BN_/\t/g' | sed 's/[.]dd_/\t/g' | sed 's/[.]a/\ta/g' | sed 's/[.]b/\tb/g' > ${summary}
	    if [[ ${noisedirname} == ${DD_DIR} ]];
	    then
		cat  ${noisedirname}/results/${method}/*/*${addon}.txt | cut -f1,2,3 | sort -k3 -n |sed 's/Depth_//g' | sed 's/[.]G/\tG/g' | sed 's/[.]K/\tK/g' | sed 's/[.]I/\tI/g' | sed 's/[.]H/\tH/g' | sed 's/[.]NH/\tNH/g' | sed 's/[.]EN_/\t/g' | sed 's/[.]NN_/\t/g' | sed 's/[.]eps_0[.]9//g' | sed 's/[.]BN_/\t/g' | sed 's/[.]dd_/\t/g' | sed 's/[.]a/\ta/g' | sed 's/[.]b/\tb/g' > ${summary}
		#cat  ${noisedirname}/${method}/*/*scores.txt |  cut -f1,2,3 | sed 's/D//g' | sed 's/[.]TM/\t/g' | sed 's/[.]S/\t/g' |sed 's/[.]EN/\t/g' | sed 's/[.]NN/\t/g'| sed 's/[.]BN/\t/g' | sed 's/[.]dd_/\t/g' | sed 's/[.]a/\ta/g' | sed 's/[.]b/\tb/g' | awk '{ddsame="different"}{if ($7==$15) ddsame="same"}{if ($3==$11) print $1"\t"$2"\t"$4"\t"$6"\t"ddsame"\t"$3"\t"$11"\t"$17}' > ${summary}
	    fi
	    
	    if [[ ${noisedirname} == ${NONREP_DIR} ]];
            then
		cat  ${noisedirname}/results/${method}/*/*${addon}.txt | cut -f1,2,3 | sort -k3 -n |sed 's/Depth_//g' | sed 's/[.]G/\tG/g' | sed 's/[.]K/\tK/g' | sed 's/[.]I/\tI/g' | sed 's/[.]H/\tH/g' | sed 's/[.]NH/\tNH/g' | sed 's/[.]EN_/\t/g' | sed 's/[.]NN_/\t/g' | sed 's/[.]eps_0[.]9//g' | sed 's/[.]BN_/\t/g' | sed 's/[.]dd_/\t/g' | sed 's/[.]a/\ta/g' | sed 's/[.]b/\tb/g' > ${summary}
		#sed 's/[.]/\t/g' | sed 's/Depth_//g' | sed 's/EN_//g' | sed 's/NN_//g' | sed 's/BN_//g' | sed 's/dd_//g' | awk '{print $1"\t"$2"\t"$13"\t"$3"."$4"\t"$7"."$8"\t"$9"\t"$10"\t"$21"\t"$23"\t"$11"."$12}' > ${summary}
		#cut -f1,2,3 | sed 's/D//g' | sed 's/[.]TM/\t/g' | sed 's/[.]S/\t/g' |sed 's/[.]EN/\t/g' | sed 's/[.]NN/\t/g'| sed 's/[.]BN/\t/g' | sed 's/[.]a/\ta/g' | sed 's/[.]b/\tb/g' | awk '{rep="nonrep"}{if ($3==$10) rep="biorep"}{print $1"\t"$2"\t"$4"\t"$6"\t"rep"\t"$3"\t"$15}' > ${summary}
	    fi
	done
    done
fi








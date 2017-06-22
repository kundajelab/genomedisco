
desired_action=$1

#=========================

CODE=/srv/gsfs0/projects/snyder/oursu/software/git/public_genomedisco/genomedisco
bashrc=${CODE}/scripts/bashrc_genomedisco
source ${bashrc}

#=========================

resolution=40000
maxdist=5000000
tadfile=/srv/gsfs0/projects/kundaje/users/oursu/3d/LA/tads/tads.chr21.merged.gz
nodefile=/ifs/scratch/oursu/encode_highres/results/res40000/data/nodes/Nodes.w40000.bed.gz.chr21.gz
realdatafile=/srv/gsfs0/projects/kundaje/users/oursu/3d/LA/merged_nodups/processed_data/results/results/ProcessedData/HIC003/HIC003.chr21.original.npz
tadmeansize=240000 #median of real data on chr21
intertadmeandistance=${resolution}
numsim=5

SIMULATION_DIR=/ifs/scratch/oursu/3d/paper/2017-05-30/simulations
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
eps=0.9


if [[ "${desired_action}" == "EdgeNoise" ]];
then

    #Edge noise
    for depth in 10000 100000 1000000 10000000;
    do
	edgenoise="0.0,0.1,0.25,0.5,0.75,0.9"
	nodenoise=0.0
	boundarynoise=0
	cmd="${mypython} ${CODEDIR}/genomedisco/simulations.py --outdir ${EDGE_DIR} --resolution ${resolution} --maxdist ${maxdist} --tadfile ${tadfile} --nodefile ${nodefile} --realdatafile ${realdatafile} --tadmeansize ${tadmeansize} --intertadmeandistance ${intertadmeandistance} --depth ${depth} --edgenoise ${edgenoise} --nodenoise ${nodenoise} --eps ${eps} --boundarynoise ${boundarynoise} --numsim  ${numsim}"
	s=${EDGE_DIR}/script_depth${depth}.edgenoise$(echo ${edgenoise} | sed 's/,/_/g').sh
	echo "source ${bashrc}" > ${s}
	echo $cmd >> ${s}
	chmod 755 ${s}
	echo $s
	qsub -l hostname=scg3* -o ${s}.o -e ${s}.e ${s} 
    done
fi

if [[ "${desired_action}" == "NodeNoise" ]];
then
    #Node noise
    for depth in 10000 100000 1000000 10000000; 
    do
	edgenoise=0.0
	nodenoise="0.0,0.1,0.25,0.5,0.75,0.9"
	boundarynoise=0
	cmd="${mypython} ${CODE}/genomedisco/simulations.py --outdir ${NODE_DIR} --resolution ${resolution} --maxdist ${maxdist} --tadfile ${tadfile} --nodefile ${nodefile} --realdatafile ${realdatafile} --tadmeansize ${tadmeansize} --intertadmeandistance ${intertadmeandistance} --depth ${depth} --edgenoise ${edgenoise} --nodenoise ${nodenoise} --eps ${eps} --boundarynoise ${boundarynoise} --numsim  ${numsim}"
	s=${NODE_DIR}/script_depth${depth}.nodenoise$(echo ${nodenoise} | sed 's/,/_/g').sh
	echo "source ${bashrc}" > ${s}
	echo $cmd >> ${s}
	chmod 755 ${s}
	echo $s
	qsub -l hostname=scg3* -o ${s}.o -e ${s}.e ${s}                                                                                                                 
    done
fi

if [[ "${desired_action}" == "BoundaryNoise" ]];
then
    #Boundary noise
    for depth in 10000 100000 1000000 10000000;
    do
	edgenoise=0.0
	nodenoise=0.0
	boundarynoise="0,40000,80000,120000,160000,200000"
	cmd="${mypython} ${CODE}/genomedisco/simulations.py --outdir ${B_DIR} --resolution ${resolution} --maxdist ${maxdist} --tadfile ${tadfile} --nodefile ${nodefile} --realdatafile ${realdatafile} --tadmeansize ${tadmeansize} --intertadmeandistance ${intertadmeandistance} --depth ${depth} --edgenoise ${edgenoise} --nodenoise ${nodenoise} --eps ${eps} --boundarynoise ${boundarynoise} --numsim  ${numsim}"
	s=${B_DIR}/script_depth${depth}.boundarynoise$(echo ${boundarynoise} | sed 's/,/_/g').sh
	echo "source ${bashrc}" > ${s}
	echo $cmd >> ${s}
	chmod 755 ${s}
	echo $s
	qsub -l hostname=scg3* -o ${s}.o -e ${s}.e ${s}
    done
fi

if [[ "${desired_action}" == "DistDep" ]];
then
    #diff distance dependence
    d1=/srv/gsfs0/projects/kundaje/users/oursu/3d/LA/merged_nodups/processed_data/results/results/ProcessedData/HIC001/HIC001.chr21.original.npz
    d2=/srv/gsfs0/projects/kundaje/users/oursu/3d/LA/merged_nodups/processed_data/results/results/ProcessedData/HIC025/HIC025.chr21.original.npz
    for depth in 10000 100000 1000000 10000000;
    do
	edgenoise=0.0
	nodenoise=0.0
	boundarynoise=0
	cmd="${mypython} ${CODE}/genomedisco/simulations.py --outdir ${DD_DIR} --resolution ${resolution} --maxdist ${maxdist} --tadfile ${tadfile} --nodefile ${nodefile} --realdatafile ${realdatafile} --tadmeansize ${tadmeansize} --intertadmeandistance ${intertadmeandistance} --depth ${depth} --edgenoise ${edgenoise} --nodenoise ${nodenoise} --eps ${eps} --boundarynoise ${boundarynoise} --numsim  ${numsim} --distance_dep_data ${d1},${d2} --distance_dep_tads ${tadfile},${tadfile}"
	s=${DD_DIR}/script_depth${depth}.ddHIC${HICNUM}.sh
	echo "source ${bashrc}" > ${s}
	echo $cmd >> ${s}
	chmod 755 ${s}
	echo $s
	qsub -l hostname=scg3* -o ${s}.o -e ${s}.e ${s}
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
            for sim in {0..4};
            do
		dataname_short=D${depth}.TM${tadmeansize}.S${sim}.EN${edgenoise}.NN${nodenoise}.BN${boundarynoise}
		dataname_basic=res_${resolution}.Depth_${depth}.MaxDist_${maxdist}.simulatedTADs_mean${tadmeansize}.S_${sim}.EN_${edgenoise}_eps_${eps}.NN_${nodenoise}.BN_${boundarynoise}
		echo "${dataname_short}.adelim${dataname_short}.bdelimchr21" | sed 's/delim/\t/g' >> ${edgenoise_metadata_pairs}
		for ab in a b;
		do
                    dataname=${dataname_basic}.${ab}
                    chromo="simulated"
                    echo "${dataname_short}.${ab}delim${chromo}delim${nodefile}delim${EDGE_DIR}/${dataname}.gzdelimNA" | sed 's/delim/\t/g' >> ${edgenoise_metadata}
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
            for sim in {0..4};
            do
		dataname_short=D${depth}.TM${tadmeansize}.S${sim}.EN${edgenoise}.NN${nodenoise}.BN${boundarynoise}
		dataname_basic=res_${resolution}.Depth_${depth}.MaxDist_${maxdist}.simulatedTADs_mean${tadmeansize}.S_${sim}.EN_${edgenoise}_eps_${eps}.NN_${nodenoise}.BN_${boundarynoise}
		echo "${dataname_short}.adelim${dataname_short}.bdelimchr21" | sed 's/delim/\t/g' >> ${nodenoise_metadata_pairs}
		for ab in a b;
		do
                    dataname=${dataname_basic}.${ab}
                    chromo="simulated"
                    echo "${dataname_short}.${ab}delim${chromo}delim${nodefile}delim${NODE_DIR}/${dataname}.gzdelimNA" | sed 's/delim/\t/g' >> ${nodenoise_metadata}
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
	for boundarynoise in 0 40000 80000 120000 160000 200000;
	do
            for sim in {0..4};
            do
		dataname_short=D${depth}.TM${tadmeansize}.S${sim}.EN${edgenoise}.NN${nodenoise}.BN${boundarynoise}
		dataname_basic=res_${resolution}.Depth_${depth}.MaxDist_${maxdist}.simulatedTADs_mean${tadmeansize}.S_${sim}.EN_${edgenoise}_eps_${eps}.NN_${nodenoise}.BN_${boundarynoise}
		echo "${dataname_short}.adelim${dataname_short}.bdelimchr21" | sed 's/delim/\t/g' >> ${boundarynoise_metadata_pairs}
		for ab in a b;
		do
                    dataname=${dataname_basic}.${ab}
                    chromo="simulated"
                    echo "${dataname_short}.${ab}delim${chromo}delim${nodefile}delim${B_DIR}/${dataname}.gzdelimNA" | sed 's/delim/\t/g' >> ${boundarynoise_metadata}
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
            for sim1 in {0..4};
            do
		for  ab1 in a b;
		do 
		    dataname_short1=D${depth}.TM${tadmeansize}.S${sim1}.EN${edgenoise}.NN${nodenoise}.BN${boundarynoise}
		    dataname_basic1=res_${resolution}.Depth_${depth}.MaxDist_${maxdist}.simulatedTADs_mean${tadmeansize}.S_${sim1}.EN_${edgenoise}_eps_${eps}.NN_${nodenoise}.BN_${boundarynoise}.${ab1}
		    for sim2 in {0..4};
		    do
			for ab2 in a b;
			do
			    dataname_short2=D${depth}.TM${tadmeansize}.S${sim2}.EN${edgenoise}.NN${nodenoise}.BN${boundarynoise}
			    dataname_basic2=res_${resolution}.Depth_${depth}.MaxDist_${maxdist}.simulatedTADs_mean${tadmeansize}.S_${sim2}.EN_${edgenoise}_eps_${eps}.NN_${nodenoise}.BN_${boundarynoise}.${ab2}
			
			    echo "${dataname_short1}.${ab1}delim${dataname_short2}.${ab2}delimsimulated" | sed 's/delim/\t/g' >> ${nonrep_metadata_pairs}.many
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
	for  dd1 in 0 1;
	do
	    for sim1 in {0..4};
	    do
		for  ab1 in a b;
		do
		    dataname_short1=D${depth}.TM${tadmeansize}.S${sim1}.EN${edgenoise}.NN${nodenoise}.BN${boundarynoise}.dd_${dd1}
		    dataname_basic1=res_${resolution}.Depth_${depth}.MaxDist_${maxdist}.simulatedTADs_mean${tadmeansize}.S_${sim1}.EN_${edgenoise}_eps_${eps}.NN_${nodenoise}.BN_${boundarynoise}.${ab1}.dd_${dd1}
		    for   dd2 in 0 1;
		    do
			for sim2 in ${sim1}; #$(eval echo "{$sim1..4}");
			do                            	   
			    for ab2 in a b;
			    do
				dataname_short2=D${depth}.TM${tadmeansize}.S${sim2}.EN${edgenoise}.NN${nodenoise}.BN${boundarynoise}.dd_${dd2}
				dataname_basic2=res_${resolution}.Depth_${depth}.MaxDist_${maxdist}.simulatedTADs_mean${tadmeansize}.S_${sim2}.EN_${edgenoise}_eps_${eps}.NN_${nodenoise}.BN_${boundarynoise}.${ab2}.dd_${dd2}
		    
				echo "${dataname_short1}.${ab1}delim${dataname_short2}.${ab2}delimsimulated" | sed 's/delim/\t/g' >> ${dd_metadata_pairs}.many
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
            for sim1 in {0..4};
            do
	        for ab in a b;
		do
		    dataname_short1=D${depth}.TM${tadmeansize}.S${sim1}.EN${edgenoise}.NN${nodenoise}.BN${boundarynoise}.dd_${dd1}.${ab}
		    dataname_basic1=res_${resolution}.Depth_${depth}.MaxDist_${maxdist}.simulatedTADs_mean${tadmeansize}.S_${sim1}.EN_${edgenoise}_eps_${eps}.NN_${nodenoise}.BN_${boundarynoise}.${ab}.dd_${dd1}
		    chromo="simulated"
		    echo "${dataname_short1}delim${chromo}delim${nodefile}delim${DD_DIR}/${dataname_basic1}.gzdelimNA" | sed 's/delim/\t/g' >> ${dd_metadata}
		done
	    done
	done
    done
fi

if [[ "${desired_action}" == "automatic" ]];
then
    for spot in DistanceDependence BoundaryNoise RepNonrep NodeNoise EdgeNoise;
    do
	chromosome=21
	comparisons=${SIMULATION_DIR}/${spot}/Metadata.pairs 
	samples=${SIMULATION_DIR}/${spot}/Metadata.samples
	outdir=${SIMULATION_DIR}/${spot}
	prefix=${spot}
	h=10
	step=simple
	resolution=40000
	bashrc=${CODEDIR}/scripts/bashrc_genomedisco
	nodes=${nodefile}
	mkdir -p ${outdir}
	script=${CODEDIR}/paper_analysis/2017-02-05/runMethods.sh

	${script} ${comparisons} ${outdir} ${prefix} ${h} ${samples} ${step} ${bashrc} ${nodes} ${resolution} ${chromosome} no
    done
fi


if [[ "${desired_action}" == "scores" ]];
then

    for noisedirname in ${EDGE_DIR} ${NODE_DIR} ${B_DIR} ${DD_DIR} ${NONREP_DIR};
    do
	for method in disco hicrep hic-spector;
	do
	    summary=${noisedirname}/${method}.results.txt
	    rm ${summary} 
	    cat  ${noisedirname}/${method}/*/*scores.txt |  cut -f1,2,3 | sed 's/D//g' | sed 's/[.]TM/\t/g' | sed 's/[.]S/\t/g' |sed 's/[.]EN/\t/g' | sed 's/[.]NN/\t/g'| sed 's/[.]BN/\t/g' | sed 's/[.]a//g' | cut -f1-6,13 > ${summary}
	    echo "====="
	    echo ${summary}
	    cat ${summary}
	    if [[ ${noisedirname} == ${DD_DIR} ]];
	    then
		cat  ${noisedirname}/${method}/*/*scores.txt |  cut -f1,2,3 | sed 's/D//g' | sed 's/[.]TM/\t/g' | sed 's/[.]S/\t/g' |sed 's/[.]EN/\t/g' | sed 's/[.]NN/\t/g'| sed 's/[.]BN/\t/g' | sed 's/[.]dd_/\t/g' | sed 's/[.]a/\ta/g' | sed 's/[.]b/\tb/g' | awk '{ddsame="different"}{if ($7==$15) ddsame="same"}{if ($3==$11) print $1"\t"$2"\t"$4"\t"$6"\t"ddsame"\t"$3"\t"$11"\t"$17}' > ${summary}
	    fi
	    
	    if [[ ${noisedirname} == ${NONREP_DIR} ]];
            then
		cat  ${noisedirname}/${method}/*/*scores.txt |  cut -f1,2,3 | sed 's/D//g' | sed 's/[.]TM/\t/g' | sed 's/[.]S/\t/g' |sed 's/[.]EN/\t/g' | sed 's/[.]NN/\t/g'| sed 's/[.]BN/\t/g' | sed 's/[.]a/\ta/g' | sed 's/[.]b/\tb/g' | awk '{rep="nonrep"}{if ($3==$10) rep="biorep"}{print $1"\t"$2"\t"$4"\t"$6"\t"rep"\t"$3"\t"$15}' > ${summary}

	    fi
	done
    done
fi








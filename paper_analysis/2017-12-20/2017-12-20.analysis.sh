step=$1
substep=$2

#=======================================
DATA=/ifs/scratch/oursu/paper_2017-12-20
resolutions=10000,50000,500000
MYCODE=/srv/gsfs0/projects/kundaje/users/oursu/code
MAPQ=30
#=======================================


if [[ ${step} == "setup_folder" ]];
then
    mkdir -p ${DATA}
fi

if [[ ${step} == "install_code" ]];
then
    #3DChromatin_ReplicateQC
    cd ${MYCODE}
    git clone http://github.com/kundajelab/3DChromatin_ReplicateQC
    cd 3DChromatin_ReplicateQC
    install_scripts/install_3DChromatin_ReplicateQC.sh --pathtopython /srv/gsfs0/projects/kundaje/users/oursu/code/anaconda2/mypython/bin/python --pathtor R --rlib /srv/gsfs0/projects/kundaje/users/oursu/code/Rlibs --pathtobedtools bedtools --modules R/3.4.1,bedtools/2.26.0
    #for hic file extraction
    mkdir -p ${MYCODE}/juicer_tools
    cd ${MYCODE}/juicer_tools
    wget http://hicfiles.tc4ga.com.s3.amazonaws.com/public/juicer/juicer_tools.1.7.5_linux_x64_jcuda.0.8.jar
    #genome_utils
    git clone https://github.com/oursu/genome_utils
    
fi

if [[ ${step} == "get_Rao_data" ]];
then
    if [[ ${substep} == "hic_files" ]];
    then
	mkdir -p ${DATA}/Rao_data/hic
	mv /ifs/scratch/oursu/data/hic/* ${DATA}/Rao_data/hic/
	for f in $(ls ${DATA}/Rao_data/hic/*gz);
	do
	    gunzip $f
	done
    fi
    if [[ ${substep} == "individual_reads" ]];
    then
	mkdir -p ${DATA}/Rao_data/counts
	mv /srv/gsfs0/projects/kundaje/users/oursu/3d/LA/merged_nodups/data/* ${DATA}/Rao_data/counts/
    fi
fi

if [[ ${step} == "Rao" ]];
then
    if [[ ${substep} == "resolutions_hic_files" ]];
    then
	for f in $(ls ${DATA}/Rao_data/hic/*hic);
	do
	    s=${DATA}/Rao_data/hic/split_$(basename ${f}).sh
	    echo "module load java/latest" > ${s}
	    for res in $(echo ${resolutions} | sed 's/,/ /g');
	    do
		mkdir -p ${DATA}/Rao_data/hic/res${res}
		for chromo in $(echo "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X" | sed 's/,/ /g');
		do
		    out=${DATA}/Rao_data/hic/res${res}/$(basename ${f}).res${res}.chr${chromo}.gz
		    echo "/usr/java/latest/bin/java -jar ${MYCODE}/juicer_tools/juicer_tools.1.7.5_linux_x64_jcuda.0.8.jar dump observed NONE ${f} ${chromo} ${chromo} BP ${res} ${out}.f" >> ${s}
		    echo "zcat -f ${out}.f | awk -v chromosome=${chromo} '{print chromosome\"\t\""'$'"1\"\t\"chromosome\"\t\""'$'"2\"\t\""'$'"3}' | gzip > ${out}" >> ${s}
		    echo "rm ${out}.f" >> ${s}
		done
	    done
            chmod 755 ${s}
	    $s
	done
	#do the gm separately
	for chromo in $(echo "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X" | sed 's/,/ /g');
	do
	    echo ${chromo}
	    zcat -f ${DATA}/Rao_data/hic/GM12878_combined/10kb_resolution_intrachromosomal/chr${chromo}/MAPQGE30/chr${chromo}_10kb.RAWobserved | awk -v chromosome=${chromo} '{print chromosome"\t"$1"\t"chromosome"\t"$2"\t"$3}' | gzip > ${DATA}/Rao_data/hic/res10000/GM12878_combined.res10000.chr${chromo}.gz
	done
	for chromo in $(echo "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X" | sed 's/,/ /g');
        do
            echo ${chromo}
            zcat -f ${DATA}/Rao_data/hic/GM12878_combined/50kb_resolution_intrachromosomal/chr${chromo}/MAPQGE30/chr${chromo}_50kb.RAWobserved | awk -v chromosome=${chromo} '{print chromosome"\t"$1"\t"chromosome"\t"$2"\t"$3}' | gzip > ${DATA}/Rao_data/hic/res50000/GM12878_combined.res50000.chr${chromo}.gz
        done
	for chromo in $(echo "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X" | sed 's/,/ /g');
        do
            echo ${chromo}
            zcat -f ${DATA}/Rao_data/hic/GM12878_combined/500kb_resolution_intrachromosomal/chr${chromo}/MAPQGE30/chr${chromo}_500kb.RAWobserved | awk -v chromosome=${chromo} '{print chromosome"\t"$1"\t"chromosome"\t"$2"\t"$3}' | gzip > ${DATA}/Rao_data/hic/res500000/GM12878_combined.res500000.chr${chromo}.gz
        done

	for cellline in GM12878_combined HMEC HUVEC IMR90 K562 KBM7 NHEK;
	do
	    echo ${cellline}
	    for res in $(echo ${resolutions} | sed 's/,/ /g');
	    do
		zcat -f ${DATA}/Rao_data/hic/res${res}/*${cellline}* | gzip > ${DATA}/Rao_data/hic/res${res}/${cellline}.res${res}.gz
	    done
	done
    fi

    if [[ ${substep} == "resolutions_reads_files" ]];
    then
	for f in $(ls ${DATA}/Rao_data/counts/*gz);
	do
	    for res in $(echo ${resolutions} | sed 's/,/ /g');
            do 
		mkdir -p ${DATA}/Rao_data/counts/res${res}
		#once for the whole genome
		out=${DATA}/Rao_data/counts/res${res}/$(basename ${f} | sed 's/_merged_nodups[.]txt[.]gz[.]//g' | sed 's/_/\t/g' | cut -f2).res${res}.gz
		s=${out}.sh
		echo ${out}
		echo "${MYCODE}/genome_utils/3Dutils/LA_reads_to_n1n2value_bins.sh ${f} ${out} ${MAPQ} intra ${res}" > ${s}
		echo "rm ${s}" >> ${s}
		chmod 755 ${s}
		qsub -o ${s}.o -e ${s}.e ${s}
	    done
	done
    fi
fi

metadata=$1
outdir=$2
prefix=$3
h=$4
node_metadata=$5
desired_action=$6
bashrc=$7
nodefile=$8
resolution=$9
chromosome=${10}
names=${11}
remove=${12}
CODE=/srv/gsfs0/projects/snyder/oursu/software/git/public_genomedisco/genomedisco

#metadata=/srv/gsfs0/projects/snyder/oursu/3d/simulations/EdgeNoise/Metadata.pairs
#outdir=/srv/gsfs0/projects/snyder/oursu/software/git/public_genomedisco/genomedisco/paper_analysis/2017-02-05/Testing
#prefix=prefix
#h=5
#node_metadata=/srv/gsfs0/projects/snyder/oursu/3d/simulations/EdgeNoise/Metadata.samples
#desired_action=simple
#bashrc=/srv/gsfs0/projects/snyder/oursu/software/git/public_genomedisco/genomedisco/scripts/bashrc_genomedisco
#nodefile=/srv/gsfs0/projects/kundaje/users/oursu/3d/encode_meeting1/AllChrAnon/processed/nodes/Nodes.w40000.chr21.bed.gz
#resolution=40000
#chromosome=21

if [[ "${desired_action}" == "simple" ]];
then
    source ${bashrc}
    echo ${nodefile} 
    zcat -f ${nodefile} | awk -v res=${resolution} '{e=$2+res/2}{print $1"\t"$2"\t"$3"\t"e}' | head
    zcat -f ${nodefile} | awk -v res=${resolution} '{e=$2+res/2}{print $1"\t"$2"\t"$3"\t"e}' | gzip > ${outdir}/nodes.res${resolution}.${chromosome}.plusResOver2.gz
    while read line
    do
	read -a items <<< "$line"
	firstItem=${items[0]}
	secondItem=${items[1]}
	pair=${firstItem}_${secondItem}
	f1=$(zcat -f ${node_metadata} | cut -f1-4 | grep ${firstItem} | cut -f4)
	f2=$(zcat -f ${node_metadata} | cut -f1-4 |  grep ${secondItem} | cut -f4)
	
	mkdir -p ${outdir}/hicrep/${pair}
	outname=${outdir}/hicrep/${pair}/chr${chromosome}.${firstItem}.vs.${secondItem}.txt
	s=${outname}.sh
	hicrepcode="${CODE}/paper_analysis/2017-02-05/other_methods/2017-01-12_run_hicrep.R"
	cmd="Rscript ${hicrepcode} ${f1} ${f2} ${outname} 1 2 3 5000000 ${resolution} ${nodefile} ${names} ${h} ${firstItem} ${secondItem}"
	echo "source ${bashrc}" > ${s}
	echo ${cmd} >> ${s}
	chmod 755 ${s}
	qsub -l h_vmem=20G -l hostname=scg3* -o ${s}.o -e ${s}.e ${s}
    
	mkdir -p ${outdir}/hic-spector/${pair}
	outname=${outdir}/hic-spector/${pair}/chr${chromosome}.${firstItem}.vs.${secondItem}.txt
	s=${outname}.sh
	hicspectorcode="${CODE}/paper_analysis/2017-02-05/other_methods/run_hicspector.jl"
	cmd="julia ${hicspectorcode} ${outname}.f1 ${outname}.f2 ${outname} ${resolution} ${chromosome} ${firstItem} ${secondItem}"
	echo "source ${bashrc}" > ${s}
	echo "zcat -f ${f1} | sed 's/[.]0//g' | awk '{print "'$'"1\"\t\""'$'"2\"\t\""'$'"3\"\n\""'$'"2\"\t\""'$'"1\"\t\""'$'"3}' | sort | uniq  > ${outname}.f1" >> ${s}
	echo "zcat -f ${f2} | sed 's/[.]0//g' | awk '{print "'$'"1\"\t\""'$'"2\"\t\""'$'"3\"\n\""'$'"2\"\t\""'$'"1\"\t\""'$'"3}' | sort | uniq > ${outname}.f2" >> ${s}
	echo ${cmd} >> ${s}
	echo "rm ${outname}.f1 ${outname}.f2" >> ${s}
	chmod 755 ${s}
	qsub -l h_vmem=20G -l hostname=scg3* -o ${s}.o -e ${s}.e ${s}
	
	mkdir -p ${outdir}/disco/${pair}
        echo  ${outdir}/disco
	outname=${outdir}/disco/${pair}/chr${chromosome}.${firstItem}.vs.${secondItem}.txt
        s=${outname}.sh
	if [[ ${remove} == "remove" ]];
	then
	    f1new=${outname}.1.gz
	    f2new=${outname}.2.gz
	    #zcat -f ${f1} | cut -f2,4,5  | gzip > ${f1new}
	    #zcat -f ${f2} | cut -f2,4,5  | gzip > ${f2new}
	    f1=${f1new}
	    f2=${f2new}
	fi
	echo ${mypython}
	disco_nodes=${outdir}/nodes.res${resolution}.${chromosome}.plusResOver2.gz
	if [[ ${names} == "use" ]];
	then
	    disco_nodes=${nodefile}
	fi
	echo "source ${bashrc}" > ${s}
	echo "${mypython} ${CODE}/genomedisco/compute_reproducibility.py --m1 ${f1} --m2 ${f2} --m1name ${firstItem} --m2name ${secondItem} --node_file ${disco_nodes} --outdir ${outdir}/disco/${pair}/ --outpref chr${chromosome} --m_subsample NA --approximation 40000 --norm sqrtvc --method RandomWalks --tmin 3 --tmax 3 --concise_analysis" >> ${s}
	chmod 755 ${s}
	#echo ${s}
	qsub -o ${s}.o -e ${s}.e ${s}
	echo "========"
    done < ${metadata}
fi

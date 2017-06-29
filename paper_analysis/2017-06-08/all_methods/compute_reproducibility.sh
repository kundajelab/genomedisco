 
usage()
{
cat "usage: `basename $0` options
Contact: Oana Ursu oursu@stanford.edu
Script for computing reproducibility using genomedisco,hicrep,hic-spector

OPTIONS:
   -h     Show this message and exit
   -o name of the output directory
   -n file containing the genomic bins used. The format is chromosome, start, end, name. The name of the bin must match the entries in the contact maps.
   -s file listing all the samples. The format is sample name, file that contains the contacts
   -p file listing all pairs. The format is sample name 1, sample name 2, chromosome
   -a action to be performed. Possible options are: compute, summarize.
   -b bashrc file. Must contain the following variables:hicrepcode ( location of the wrapper around hicrep),hic-spectorcode, genomediscocode 
   -r resolution, in base pairs
   -m metadata. This file lists any parameters associated with running these methods. The format of the file is method, parameter, parameter value, tab separated. The parameters for genomedisco are tmin(minimum steps of random walk),tmax, normalization (can be sqrtvc). The parameters for hicrep are h(smoothing size), maxdist(maximum distance to consider). The parameters for hic-spector are n(the number of eigenvectors to use).
   -j  whether to run these interactively or submit as jobs, using sge.  To submit jobs, set this to 'sge',  otherwise set to 'not_sge' (DEFAULT)
   -c chromosome
"
} 

memory="-l h_vmem=20G"
#memory=""
outdir='outdir'
sge='not_parallel'
while getopts "ho:n:s:p:a:b:r:m:j:c:" opt
do
    case $opt in
	h)
	        usage; exit;;
	o)
	        outdir=$OPTARG;;
	n)
	        bins=$OPTARG;;
	s)
	        samples=$OPTARG;;
	p)
	        pairs=$OPTARG;;
	a)
	        action=$OPTARG;;
	b)
	        bashrc=$OPTARG;;
	r)
	        resolution=$OPTARG;;
	m)
	        metadata=$OPTARG;;
	j)
	        sge=$OPTARG;;
	c)
	        chromosome=$OPTARG;;
	?)
	        usage
		    exit 1;;
    esac    
done

if [[ "${action}" == "compute" ]];
then
    source ${bashrc}
    #zcat -f ${bins} | awk -v res=${resolution} '{mid=$2+res/2}{print $1"\t"$2"\t"$3"\t"mid}' | gzip > ${outdir}/$(basename ${bins}).mid.gz

    while read line
    do
	read -a items <<< "$line"
	firstItem=${items[0]}
	secondItem=${items[1]}
	#chromosome=${items[2]}
	chromosome_number=$(echo ${chromosome} | sed 's/chr//g')
	pair=${firstItem}_${secondItem}
	f1=$(zcat -f ${samples} | cut -f1-2 | grep ${firstItem} | cut -f2)
	f2=$(zcat -f ${samples} | cut -f1-2 |  grep ${secondItem} | cut -f2)
	
	mkdir -p ${outdir}/results/hicrep/${pair}
	outname=${outdir}/results/hicrep/${pair}/hicrep.${chromosome}.${firstItem}.vs.${secondItem}.txt
	s=${outname}.sh
	h=$(zcat -f ${metadata} | grep "hicrep"  | grep -w h | sed 's/_/\t/g' | cut -f3)
	maxdist=$(zcat -f ${metadata} | grep "hicrep"  | grep -w maxdist | sed 's/_/\t/g' | cut -f3)
	
	cmd="Rscript ${hicrepcode} ${f1} ${f2} ${outname} ${maxdist} ${resolution} ${bins} ${h} ${firstItem} ${secondItem}"
	echo "source ${bashrc}" > ${s}
	echo ${cmd} >> ${s}
	chmod 755 ${s}
	#if [[ ${sge} == "sge" ]];
	#then
	#    qsub ${memory} -o ${s}.o -e ${s}.e ${s}
	#else
	#    ${s}
	#fi

	mkdir -p ${outdir}/results/hic-spector/${pair}
	outname=${outdir}/results/hic-spector/${pair}/hic-spector.${chromosome}.${firstItem}.vs.${secondItem}.txt
	n=$(zcat -f ${metadata} | grep "hic-spector"  | grep -w n | sed 's/_/\t/g' | cut -f3)
	s=${outname}.sh
	cmd="julia ${hicspectorcode} ${outname}.f1 ${outname}.f2 ${outname} ${resolution} ${chromosome_number} ${firstItem} ${secondItem} ${n} ${spector}/HiC_spector.jl"
	echo "source ${bashrc}" > ${s}
	echo "zcat -f ${f1} | sed 's/[.]0//g' | awk '{print "'$'"1\"\t\""'$'"2\"\t\""'$'"3\"\n\""'$'"2\"\t\""'$'"1\"\t\""'$'"3}' | sort | uniq  > ${outname}.f1" >> ${s}
	echo "zcat -f ${f2} | sed 's/[.]0//g' | awk '{print "'$'"1\"\t\""'$'"2\"\t\""'$'"3\"\n\""'$'"2\"\t\""'$'"1\"\t\""'$'"3}' | sort | uniq > ${outname}.f2" >> ${s}
	echo ${cmd} >> ${s}
	echo "rm ${outname}.f1 ${outname}.f2" >> ${s}
	chmod 755 ${s}
	if [[ ${sge} == "sge" ]];
        then
	    qsub ${memory} -o ${s}.o -e ${s}.e ${s}
	else
            ${s}
        fi
	
	mkdir -p ${outdir}/results/genomedisco/${pair}
	tmin=$(zcat -f ${metadata} | grep "genomedisco"  | grep -w tmin | sed 's/_/\t/g' | cut -f3)
	tmax=$(zcat -f ${metadata} | grep "genomedisco"  | grep -w tmax | sed 's/_/\t/g' | cut -f3)
	normalization=$(zcat -f ${metadata} | grep "genomedisco"  | grep -w normalization | sed 's/_/\t/g' | cut -f3)
	outname=${outdir}/results/genomedisco/${pair}/genomedisco.${chromosome}.${firstItem}.vs.${secondItem}.txt
        s=${outname}.sh
	echo "source ${bashrc}" > ${s}
	echo "${mypython} ${genomedisco}/genomedisco/compute_reproducibility.py --m1 ${f1} --m2 ${f2} --m1name ${firstItem} --m2name ${secondItem} --node_file ${bins} --outdir ${outdir}/results/genomedisco/${pair}/ --outpref ${chromosome} --m_subsample NA --approximation 40000 --norm ${normalization} --method RandomWalks --tmin ${tmin} --tmax ${tmax} --concise_analysis" >> ${s}
	chmod 755 ${s}
	#if [[ ${sge} == "sge" ]];
        #then
        #    qsub -l h_vmem=50G -o ${s}.o -e ${s}.e ${s}
	#else
        #    ${s}
        #fi
    done < ${pairs}
fi


if [[ "${action}" == "scorelist" ]];
then
    source ${bashrc}
    mkdir -p ${outdir}/scorelist
    chromosomes=$(zcat -f ${pairs} | cut -f3 | sort | uniq) 
    for chromosome in ${chromosomes};
    do
	scores=${outdir}/scorelist/${prefix}.${chromosome}.txt 
	echo "#m1_m2_chromosome_genomedisco_hicrep_hic-spector" | sed 's/_/\t/g' > ${scores}
    done
    while read line
    do
        read -a items <<< "$line"
        firstItem=${items[0]}
        secondItem=${items[1]}
        chromosome=${items[2]}
        chromosome_number=$(echo ${chromosome} | sed 's/chr//g')
	pair=${firstItem}_${secondItem}
	scores=${outdir}/scorelist/${prefix}.${chromosome}.txt 
	
	d=$(cat ${outdir}/genomedisco/${pair}/genomedisco.${chromosome}.${firstItem}.vs.${secondItem}.txt | cut -f3)
	r=$(cat ${outdir}/hicrep/${pair}/hicrep.${chromosome}.${firstItem}.vs.${secondItem}.txt | cut -f3)
	s=$(cat ${outdir}/hic-spector/${pair}/hic-spector.${chromosome}.${firstItem}.vs.${secondItem}.txt | cut -f3)
	
	echo "${firstItem}_${secondItem}_${chromosome}_${d}_${r}_${s}" | sed 's/_/\t/g' >> ${scores}
    done < ${pairs}
fi
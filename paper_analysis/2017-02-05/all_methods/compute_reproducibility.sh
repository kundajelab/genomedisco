 
usage()
{
cat "usage: `basename $0` options
Contact: Oana Ursu oursu@stanford.edu
Script for computing reproducibility using genomedisco,hicrep,hic-spector

OPTIONS:
   -h     Show this message and exit
   -o name of the output directory
   -p prefix for output files
   -n file containing the genomic bins used. The format is chromosome, start, end, name. The name of the bin must match the entries in the contact maps.
   -s file listing all the samples. The format is sample name, file that contains the contacts
   -c file listing all comparisons. The format is sample name 1, sample name 2, chromosome
   -a action to be performed. Possible options are: compute, summarize.
   -b bashrc file. Must contain the following variables:hicrepcode ( location of the wrapper around hicrep),hic-spectorcode, genomediscocode 
   -r resolution, in base pairs
   -m metadata. This file lists any parameters associated with running these methods. The format of the file is method, parameter, parameter value, tab separated. The parameters for genomedisco are tmin(minimum steps of random walk),tmax, normalization (can be sqrtvc). The parameters for hicrep are h(smoothing size), maxdist(maximum distance to consider). The parameters for hic-spector are n(the number of eigenvectors to use).
"
} 

outdir='outdir'
prefix='reproducibility'
while getopts "ho:p:n:s:c:a:b:r:m:" opt
do
    case $opt in
	h)
	        usage; exit;;
	o)
	        outdir=$OPTARG;;
	p)
	        prefix=$OPTARG;;
	n)
	        bins=$OPTARG;;
	s)
	        samples=$OPTARG;;
	c)
	        pairs=$OPTARG;;
	a)
	        action=$OPTARG;;
	b)
	        bashrc=$OPTARG;;
	r)
	        resolution=$OPTARG;;
	m)
	        metadata=$OPTARG;;
	?)
	        usage
		    exit 1;;
    esac    
done

if [[ "${action}" == "compute" ]];
then
    source ${bashrc}
    zcat -f ${bins} | awk -v res=${resolution} '{mid=$2+res/2}{print $1"\t"$2"\t"$3"\t"mid}' | gzip > ${outdir}/$(basename ${bins}).mid.gz

    while read line
    do
	read -a items <<< "$line"
	firstItem=${items[0]}
	secondItem=${items[1]}
	chromosome=${items[2]}
	chromosome_number=$(echo ${chromosome} | sed 's/chr//g')
	pair=${firstItem}_${secondItem}
	f1=$(zcat -f ${samples} | awk -v chrom=${chromosome} '{if ($3==chrom) print $0}' | cut -f1-2 | grep ${firstItem} | cut -f2)
	f2=$(zcat -f ${samples} | awk -v chrom=${chromosome} '{if ($3==chrom) print $0}' | cut -f1-2 |  grep ${secondItem} | cut -f2)
	
	mkdir -p ${outdir}/hicrep/${pair}
	outname=${outdir}/hicrep/${pair}/hicrep.${chromosome}.${firstItem}.vs.${secondItem}.txt
	s=${outname}.sh
	h=$(zcat -f ${metadata} | grep "hicrep"  | grep -w h | sed 's/_/\t/g' | cut -f3)
	maxdist=$(zcat -f ${metadata} | grep "hicrep"  | grep -w maxdist | sed 's/_/\t/g' | cut -f3)
	
	cmd="Rscript ${hicrepcode} ${f1} ${f2} ${outname} ${maxdist} ${resolution} ${bins} ${h} ${firstItem} ${secondItem}"
	echo "source ${bashrc}" > ${s}
	echo ${cmd} >> ${s}
	chmod 755 ${s}
	#qsub -l h_vmem=20G -l hostname=scg3* -o ${s}.o -e ${s}.e ${s}
    
	mkdir -p ${outdir}/hic-spector/${pair}
	outname=${outdir}/hic-spector/${pair}/hic-spector.${chromosome}.${firstItem}.vs.${secondItem}.txt
	n=$(zcat -f ${metadata} | grep "hic-spector"  | grep -w n | sed 's/_/\t/g' | cut -f3)

	s=${outname}.sh
	cmd="julia ${hicspectorcode} ${outname}.f1 ${outname}.f2 ${outname} ${resolution} ${chromosome_number} ${firstItem} ${secondItem} ${n} ${spector}/HiC_spector.jl"
	echo "source ${bashrc}" > ${s}
	echo "zcat -f ${f1} | sed 's/[.]0//g' | awk '{print "'$'"1\"\t\""'$'"2\"\t\""'$'"3\"\n\""'$'"2\"\t\""'$'"1\"\t\""'$'"3}' | sort | uniq  > ${outname}.f1" >> ${s}
	echo "zcat -f ${f2} | sed 's/[.]0//g' | awk '{print "'$'"1\"\t\""'$'"2\"\t\""'$'"3\"\n\""'$'"2\"\t\""'$'"1\"\t\""'$'"3}' | sort | uniq > ${outname}.f2" >> ${s}
	echo ${cmd} >> ${s}
	echo "rm ${outname}.f1 ${outname}.f2" >> ${s}
	chmod 755 ${s}
	#qsub -l h_vmem=20G -l hostname=scg3* -o ${s}.o -e ${s}.e ${s}
	
	mkdir -p ${outdir}/genomedisco/${pair}
	tmin=$(zcat -f ${metadata} | grep "genomedisco"  | grep -w tmin | sed 's/_/\t/g' | cut -f3)
	tmax=$(zcat -f ${metadata} | grep "genomedisco"  | grep -w tmax | sed 's/_/\t/g' | cut -f3)
	normalization=$(zcat -f ${metadata} | grep "genomedisco"  | grep -w normalization | sed 's/_/\t/g' | cut -f3)
	outname=${outdir}/genomedisco/${pair}/genomedisco.${chromosome}.${firstItem}.vs.${secondItem}.txt
        s=${outname}.sh
	echo "source ${bashrc}" > ${s}
	echo "${mypython} ${genomedisco}/genomedisco/__main__.py --m1 ${f1} --m2 ${f2} --m1name ${firstItem} --m2name ${secondItem} --node_file ${bins} --outdir ${outdir}/genomedisco/${pair}/ --outpref ${chromosome} --m_subsample NA --approximation 40000 --norm ${normalization} --method RandomWalks --tmin ${tmin} --tmax ${tmax} --concise_analysis" >> ${s}
	chmod 755 ${s}
	#qsub -l h_vmem=20G -l hostname=scg3* -o ${s}.o -e ${s}.e ${s}
	echo "========"
    done < ${pairs}
fi

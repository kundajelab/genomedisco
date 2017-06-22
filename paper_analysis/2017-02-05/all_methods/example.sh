#==================
# Specify step
#==================
step=$1

#==================
# setup environment
#==================
codebase=/srv/gsfs0/projects/snyder/oursu/software/git/public_genomedisco/genomedisco
bashrc=${codebase}/paper_analysis/2017-02-05/all_methods/methods_bashrc
source ${bashrc}

#==================
# set up samples to compare
#==================
chrSizes=${codebase}/paper_analysis/2017-02-05/all_methods/chrSizes
outdir=${code}/test
mkdir -p ${outdir}
prefix=prefix
resolution=40000
samples=${code}/example_samples.txt
comparisons=${code}/example_comparisons.txt
echo "HIC001 ${code}/HIC001.res40000.gz" | tr " " "\t" > ${samples}
echo "HIC002 ${code}/HIC002.res40000.gz" | tr " " "\t" >> ${samples}
echo "HIC050 ${code}/HIC050.res40000.gz" | tr " " "\t" >> ${samples}
echo "HIC001 HIC002" | tr " " "\t" > ${comparisons}
echo "HIC001 HIC050" | tr " " "\t" >> ${comparisons}
echo "HIC002 HIC050" | tr " " "\t" >> ${comparisons}

#==================
# create genomic bins
#==================
if [[ ${step} == "bin" ]];
then
    bedtools makewindows -i winnum -w ${resolution} -s ${resolution} -g ${chrSizes} | awk '{print $1"\t"$2"\t"$3"\t"$2}' | gzip > ${outdir}/${prefix}.bins.tmp.gz
    #we'll only use chr21 and chr22
    zcat -f ${outdir}/${prefix}.bins.tmp.gz | grep chr21 > ${outdir}/${prefix}.bins.tmp.gz_chr21
    zcat -f ${outdir}/${prefix}.bins.tmp.gz | grep chr22 > ${outdir}/${prefix}.bins.tmp.gz_chr22
    zcat -f ${outdir}/${prefix}.bins.tmp.gz_chr21 ${outdir}/${prefix}.bins.tmp.gz_chr22 | gzip > ${outdir}/${prefix}.bins.gz
    rm ${outdir}/${prefix}.bins.tmp.gz_chr21 ${outdir}/${prefix}.bins.tmp.gz_chr22 ${outdir}/${prefix}.bins.tmp.gz
fi 

#==================
# split files by chromosome
#==================
if [[ ${step} == "split" ]];
then
    ${mypython} ${genomedisco}/genomedisco/__main__.py split --datatype hic --metadata_samples ${samples} --nodes ${outdir}/${prefix}.bins.gz --outdir ${outdir}
fi

#==================
# create files describing the samples and the desired comparisons.  for speed, we will do a subset of the chromosomes,chr21,chr22
#==================
if [[ ${step} == "metadata" ]];
then
    for chromosome in $(zcat -f ${outdir}/data/metadata/chromosomes.gz | sed 's/\n/ /g' | sed 's/chr//g');
    do
	zcat -f ${samples} | awk -v out=${outdir} -v chrom=${chromosome} '{print $1"\t"out"/data/edges/"$1"/"$1".chr"chrom".gz\tchr"chrom}' > ${code}/example_samples.chr${chromosome}.txt
	zcat -f ${comparisons} | awk -v chrom=${chromosome} '{print $1"\t"$2"\tchr"chrom}' > ${code}/example_comparisons.chr${chromosome}.txt
    done 
fi

if [[ ${step} == "compute" ]];
then
    for chromosome in $(zcat -f ${outdir}/data/metadata/chromosomes.gz | sed 's/\n/ /g' | sed 's/chr//g');
    do 
	echo ${chromosome}
	bins=${outdir}/data/nodes/nodes.chr${chromosome}.gz
	samples=${code}/example_samples.chr${chromosome}.txt
	comparisons=${code}/example_comparisons.txt
	parameters=${code}/example_parameters.txt
	action=compute
        ${code}/compute_reproducibility.sh -o ${outdir} -n ${bins} -s ${samples} -p ${comparisons} -a ${action} -b ${bashrc} -r ${resolution} -m ${parameters} -j not_parallel -c chr${chromosome}
    done
fi

if [[ ${step} == "scorelist" ]];
then
    for chromosome in $(zcat -f ${outdir}/data/metadata/chromosomes.gz | sed 's/\n/ /g' | sed 's/chr//g');
    do
    bins=${outdir}/data/nodes/nodes.chr${chromosome}.gz
    samples=${code}/example_samples.chr${chromosome}.txt
    comparisons=${code}/example_comparisons.txt
    parameters=${code}/example_parameters.txt
    action=scorelist
    ${code}/compute_reproducibility.sh -o ${outdir} -n ${bins} -s ${samples} -c ${comparisons} -a ${action} -b ${bashrc} -r ${resolution} -m ${parameters} -j not_parallel
    done
fi
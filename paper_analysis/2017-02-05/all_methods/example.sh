#==================
# Specify step
#==================
step=$1

#==================
# setup environment
#==================
bashrc=/srv/gsfs0/projects/snyder/oursu/software/git/public_genomedisco/genomedisco/paper_analysis/2017-02-05/all_methods/methods_bashrc
source ${bashrc}

#==================
# set up samples to compare
#==================
chrSizes=/srv/gsfs0/projects/snyder/oursu/software/git/public_genomedisco/genomedisco/paper_analysis/2017-02-05/all_methods/chrSizes
outdir=${code}/test
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
    bedtools makewindows -i winnum -w ${resolution} -s ${resolution} -g ${chrSizes} | awk '{print $1"\t"$2"\t"$3"\t"$2}' | gzip > ${outdir}/${prefix}.bins.gz
fi 

#==================
# split files by chromosome
#==================
if [[ ${step} == "split" ]];
then
    ${genomedisco}/scripts/splitByChromosome.sh -t hic -i ${samples} -n ${outdir}/${prefix}.bins.gz -j sge -o ${outdir}/${prefix}
fi

#==================
# create files describing the samples and the desired comparisons.  for speed, we will do a subset of the chromosomes,chr10,chr11
#==================
if [[ ${step} == "metadata" ]];
then
    for chromosome in 10 11;
    do
	zcat -f ${samples} | awk -v out=${outdir} -v pref=${prefix} -v chrom=${chromosome} '{print $1"\t"out"/"pref"/data/edges/"$1"/"$1".chr"chrom".gz\tchr"chrom}' > ${code}/example_samples.chr${chromosome}.txt
	zcat -f ${comparisons} | awk -v chrom=${chromosome} '{print $1"\t"$2"\tchr"chrom}' > ${code}/example_comparisons.chr${chromosome}.txt
    done 
fi

if [[ ${step} == "compute" ]];
then
    for chromosome in 10 11;
    do 
    bins=${outdir}/${prefix}/data/nodes/prefix.bins.gz.chr${chromosome}.gz
    samples=${code}/example_samples.chr${chromosome}.txt
    comparisons=${code}/example_comparisons.chr${chromosome}.txt
    parameters=${code}/example_parameters.txt
    action=compute
    ${code}/compute_reproducibility.sh -o ${outdir}/${prefix}/results -p ${prefix} -n ${bins} -s ${samples} -c ${comparisons} -a ${action} -b ${bashrc} -r ${resolution} -m ${parameters} -j sge
    done
fi

if [[ ${step} == "scorelist" ]];
then
    for chromosome in 10 11;
    do
    bins=${outdir}/${prefix}/data/nodes/prefix.bins.gz.chr${chromosome}.gz
    samples=${code}/example_samples.chr${chromosome}.txt
    comparisons=${code}/example_comparisons.chr${chromosome}.txt
    parameters=${code}/example_parameters.txt
    action=scorelist
    ${code}/compute_reproducibility.sh -o ${outdir}/${prefix}/results -p ${prefix} -n ${bins} -s ${samples} -c ${comparisons} -a ${action} -b ${bashrc} -r ${resolution} -m ${parameters} -j sge
    done
fi
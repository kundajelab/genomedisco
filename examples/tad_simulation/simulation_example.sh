mkdir simulation_output
nodefile=Nodes.w40000.chr21.bed.gz

python ../../genomedisco/simulations.py  --outdir simulation_output --resolution 40000 --maxdist 5000000 --tadfile tads.chr21.merged.gz --nodefile ${nodefile} --realdatafile HIC003.res40000.chr21.gz --tadmeansize 1000000 --intertadmeandistance 40000 --depth 1000000 --edgenoise 0.0 --nodenoise 0.0 --eps 0.9 --boundarynoise 0 --numsim 1

for f in $(ls simulation_output/*.gz | grep "a.gz\|b.gz");
do
    ls -lh $f
    python ../../genomedisco/compute_rw.py --m ${f} --node_file ${nodefile} --mname test --outdir normalized --outpref $(basename ${f}).norm_sqrtvc --tmin 1 --tmax 1 --norm sqrtvc
done
mkdir simulation_output
python ../../genomedisco/simulations.py  --outdir simulation_output --resolution 40000 --maxdist 5000000 --tadfile tads.chr21.merged.gz --nodefile Nodes.w40000.chr21.bed.gz --realdatafile HIC003.res40000.chr21.gz --tadmeansize 1000000 --intertadmeandistance 40000 --depth 1000000 --edgenoise 0.0 --nodenoise 0.0 --eps 0.9 --boundarynoise 0 --numsim 1

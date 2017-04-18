

to_include=ARGS[9];
include(to_include)
map_file1=ARGS[1];
map_file2=ARGS[2];
outname=ARGS[3]
bin_size=parse(Float64,ARGS[4]);
chr_num=parse(Int64,ARGS[5]);
M1name=ARGS[6];
M2name=ARGS[7];
hg19_info=define_hg19_genome();
r=parse(Int64,ARGS[8]); #number of eigenvectors
chr2bins,bin2loc=generate_arbitrary_mapping_files(hg19_info,bin_size);

X1=readdlm(map_file1,Int64);
X2=readdlm(map_file2,Int64);

ib=find(bin2loc[1,:].==chr_num-1);
N=length(ib);

M1=sparse(floor(Int64,X1[:,1]/bin_size)+1,floor(Int64,X1[:,2]/bin_size)+1,X1[:,3],N,N);
M2=sparse(floor(Int64,X2[:,1]/bin_size)+1,floor(Int64,X2[:,2]/bin_size)+1,X2[:,3],N,N);

M1=full(M1);
M2=full(M2);
    
evs,a1,a2=get_reproducibility(M1,M2,r);

open("$outname", "w") do f
write(f, "$M1name\t$M2name\t$evs\n")
end




import argparse
import gzip

def main():
    parser = argparse.ArgumentParser(description='')

    parser.add_argument('--nodes')
    parser.add_argument('--baits')
    parser.add_argument('--out')

    args = parser.parse_args()
    
    node_dict={}
    baits=set()
    out=gzip.open(args.out,'w')

    for line in gzip.open(args.baits,'r').readlines():
        items=line.strip().split('\t')
        chromo,start,end=items[0],items[1],items[2]
        baits.add(chromo+':'+start+'-'+end)

    for line in gzip.open(args.nodes,'r').readlines():
        items=line.strip().split('\t')
        chromo,start,end,name=items[0],items[1],items[2],items[3]
        node=chromo+':'+start+'-'+end
        if node in baits:
            out.write(chromo+'\t'+start+'\t'+end+'\t'+name+'\t'+'included'+'\n')
        else:
            out.write(chromo+'\t'+start+'\t'+end+'\t'+name+'\t'+'notIncluded'+'\n')
    out.close()

main()

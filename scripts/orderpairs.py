
import argparse

def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--file')
    parser.add_argument('--out')
    args = parser.parse_args()

    out=open(args.out,'w')
    for line in open(args.file,'r').readlines():
        items=line.strip().split()
        i1,i2=items[0],items[1]
        if i1<=i2:
            out.write(i1+'\t'+i2+'\n')
        if i1>i2:
            out.write(i2+'\t'+i1+'\n')

main()


args=commandArgs(trailingOnly=TRUE)

require("hicrep")

#=============================

f1=args[1]
f2=args[2]
out=args[3]

c1=1
c2=2
c3=3

maxdist=as.numeric(args[4])
resol=as.numeric(args[5])
nodefile=args[6]
h=as.numeric(args[7])
m1name=args[8]
m2name=args[9]

#=============================

options(scipen=999)

# read in data sets
data1=read.table(f1)
data2=read.table(f2)
d1=data1[,c(c1,c2,c3)]
d2=data2[,c(c1,c2,c3)]

# read in nodes
nodedata=read.table(nodefile)
nodedata=data.frame(nodedata,nodename=as.numeric(as.character(nodedata[,4])))
rownames(nodedata)=nodedata[,'nodename']
nodes=as.character(nodedata[,'nodename'])

# construct matrices
m1=array(0,dim=c(length(nodes),length(nodes)))
rownames(m1)=colnames(m1)=nodes
m2=array(0,dim=c(length(nodes),length(nodes)))
rownames(m2)=colnames(m2)=nodes

for (i in c(1:(dim(d1)[1]))){
    n1=d1[i,1]
    n2=d1[i,2]
    v=d1[i,3]
    m1[as.character(n1),as.character(n2)]=m1[as.character(n1),as.character(n2)]+v
    m1[as.character(n2),as.character(n1)]=m1[as.character(n2),as.character(n1)]+v
}
for (i in c(1:(dim(d2)[1]))){
    n1=d2[i,1]
    n2=d2[i,2]
    v=d2[i,3]
    m2[as.character(n1),as.character(n2)]=m2[as.character(n1),as.character(n2)]+v
    m2[as.character(n2),as.character(n1)]=m2[as.character(n2),as.character(n1)]+v
}
m1_big=data.frame(chr='chromo',n1=as.numeric(as.character(nodedata[,2])),n2=as.numeric(as.character(nodedata[,3])),m1)
m2_big=data.frame(chr='chromo',n1=as.numeric(as.character(nodedata[,2])),n2=as.numeric(as.character(nodedata[,3])),m2)

# prepare matrices
Pre_HiC <- prep(m1_big, m2_big, resol, h, maxdist)

# compute score
SCC.out = get.scc(Pre_HiC, resol, maxdist)
print(SCC.out)

# write score
scores=data.frame(M1=m1name,M2=m2name,score=SCC.out[[3]])
write.table(scores,file=out,quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')

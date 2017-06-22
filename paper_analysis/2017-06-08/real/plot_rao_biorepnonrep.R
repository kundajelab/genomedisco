
args=commandArgs(trailingOnly=TRUE)

counts_f=args[1]
scores_f=args[2]
celltypes_f=args[3]
out=args[4]

testing=function(){
	counts_f='/srv/gsfs0/projects/kundaje/users/oursu/3d/LA/merged_nodups/processed_data/TotalCountsSamplenames.gz'
	scores_f='/srv/gsfs0/projects/kundaje/users/oursu/3d/LA/merged_nodups/processed_data/UniformIsHIC053/results/results.uniform_seq_depth/GenomeWideScores/RandomWalk_v2.aggScores.forheatmap2.txt'
	celltypes_f='/srv/gsfs0/projects/kundaje/users/oursu/3d/LA/merged_nodups/processed_data/sampleCellTypes.txt'
	out='/srv/gsfs0/projects/kundaje/users/oursu/test/testplot.pdf'
	scores_f='/srv/gsfs0/projects/kundaje/users/oursu/3d/LA/merged_nodups/processed_data/UniformIsHIC053/hic-spector/hic-spector.053.results'
}

counts=read.table(counts_f)
colnames(counts)=c('Sample','Counts')
rownames(counts)=as.character(counts[,'Sample'])

celltypes=read.table(celltypes_f)
colnames(celltypes)=c('Sample','Celltype')
rownames(celltypes)=as.character(celltypes[,'Sample'])

scores=read.table(scores_f)
scores_b=data.frame(scores[,2],scores[,1],scores[,3])
colnames(scores_b)=colnames(scores)
scores=rbind(scores,scores_b)

#scatterplot by sequencing depth
scores=data.frame(scores,seqdepth=counts[as.character(scores[,1]),'Counts'],cell1=celltypes[as.character(scores[,1]),'Celltype'],
		cell2=celltypes[as.character(scores[,2]),'Celltype'])
scores=data.frame(scores,biorep=0)
scores[which(scores[,'cell1']==scores[,'cell2']),'biorep']=1
print(scores)
#roc, prc
require(PRROC)
scores_0=scores[which(as.character(scores[,'biorep'])=='0'),3]
scores_1=scores[which(as.character(scores[,'biorep'])=='1'),3]
roc=roc.curve(scores.class0=scores[,3], weights.class0=scores[,'biorep'],curve=TRUE)
pr=pr.curve(scores.class0=scores[,3], weights.class0=scores[,'biorep'],curve=TRUE)

write.table(scores,file=paste(out,'.scoresTable',sep=''),sep='\t',quote=F,row.names=F,col.names=F)

pdf(out)
require(ggplot2)
print(ggplot(scores,aes(x=seqdepth,y=V3,color=as.factor(biorep)))+xlab('Sequencing depth of reference sample')+ylab('Score')+geom_point()+theme_bw()+scale_color_manual(values = c("gray","blue")))

scores_short=scores
scores_short[,1]=gsub('HIC0','',as.character(scores_short[,1]))
scores_short[,2]=gsub('HIC0','',as.character(scores_short[,2]))
samples=unique(c(scores_short[,1],scores_short[,2]))
samples=samples[order(samples)]

m=as.matrix(array(0,dim=c(length(samples),length(samples))))
rownames(m)=colnames(m)=samples
for (i in c(1:(dim(scores_short)[1]))){
    n1=as.character(scores_short[i,1])
    n2=as.character(scores_short[i,2])
    v=scores_short[i,3]
    m[n1,n2]=v
    m[n2,n1]=v
}

require(pheatmap)
pheatmap(m,cluster_cols=FALSE,cluster_rows=FALSE,display_numbers=TRUE,number_format = "%.2f",fontsize=5)

print(roc)
print(pr)
plot(roc)
plot(pr)

dev.off()


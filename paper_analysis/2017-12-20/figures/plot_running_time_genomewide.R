

args=commandArgs(trailingOnly=TRUE)
t_file=args[1]
out=args[2]

rtimes=data.frame(read.table(t_file),timetype='real')
colnames(rtimes)=c('method','chr','m1','m2','time','timetype')
#rownames(rtimes)=paste(rtimes[,1],rtimes[,3],rtimes[,4])

#keep only ones with all chromosomes, unless it's quasar
nonquasars=which(as.character(rtimes[,'method'])!='QuASAR-Rep')
nonq=rtimes[nonquasars,]
nonq=data.frame(nonq,mim2=paste(as.character(nonq[,'m1']),as.character(nonq[,'m2'])))
pairs_with_counts=colSums(table(nonq[,c(1,7)]))
keep=names(pairs_with_counts)[which(as.numeric(as.character(pairs_with_counts))==69)]
keep=keep[1:50]

print(length(keep))
nonq=nonq[which(as.character(nonq[,7]) %in% keep),]
nonq_genome=aggregate(nonq[,'time']~nonq[,'method']+nonq[,7],FUN=sum)

#allscores=rbind(nonq_genome,quasars)
allscores=nonq_genome
colnames(allscores)=c('method','pair','time')
allscores=data.frame(allscores,minutes=as.numeric(as.character(allscores[,'time']))/60.0)

for (method in c('GenomeDISCO','HiCRep','HiC-Spector')){
    print(method)
    print(summary(allscores[which(as.character(allscores[,'method'])==method),'time']/60))
}
require(ggplot2)
pdf(out)
print(ggplot(allscores,aes(x=method,y=minutes))+geom_boxplot()+geom_jitter(alpha=0.5)
+theme_bw())
dev.off()


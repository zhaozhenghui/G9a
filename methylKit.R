library(methylKit)
file.list=list("T1D-26-1_CpG.txt","T1D-26-2_CpG.txt","T1D-26-3_CpG.txt","T1D-26-CK-1_CpG.txt","T1D-26-CK-2_CpG.txt","T1D-26-CK-3_CpG.txt","T1D-GC-1_CpG.txt","T1D-GC-2_CpG.txt","T1D-GC-3_CpG.txt","T1D-GC-CK-1_CpG.txt","T1D-GC-CK-2_CpG.txt","T1D-GC-CK-3_CpG.txt")
myobj=read(file.list,sample.id=list("T1D-26-1","T1D-26-2","T1D-26-3","T1D-26-CK-1","T1D-26-CK-2","T1D-26-CK-3","T1D-GC-1","T1D-GC-2","T1D-GC-3","T1D-GC-CK-1","T1D-GC-CK-2","T1D-GC-CK-3"),assembly="mm10",treatment=c(0,0,0,0,0,0,1,1,1,1,1,1))
meth=methylKit::unite(myobj,destrand=FALSE)

file.list=list("T1D-26-1_CpG.txt","T1D-26-2_CpG.txt","T1D-26-CK-1_CpG.txt","T1D-26-CK-2_CpG.txt","T1D-26-CK-3_CpG.txt")
myobj=read(file.list,sample.id=list("T1D-26-1","T1D-26-2","T1D-26-CK-1","T1D-26-CK-2","T1D-26-CK-3"),assembly="mm10",treatment=c(0,0,1,1,1))
meth=methylKit::unite(myobj,destrand=FALSE)

file.list=list("T1D-GC-1_CpG.txt","T1D-GC-2_CpG.txt","T1D-GC-3_CpG.txt","T1D-GC-CK-1_CpG.txt","T1D-GC-CK-2_CpG.txt","T1D-GC-CK-3_CpG.txt")
myobj=read(file.list,sample.id=list("T1D-GC-1","T1D-GC-2","T1D-GC-3","T1D-GC-CK-1","T1D-GC-CK-2","T1D-GC-CK-3"),assembly="mm10",treatment=c(0,0,0,1,1,1))
meth=methylKit::unite(myobj,destrand=FALSE)

png("correlation.png",width=1600,height=1600)
getCorrelation(meth,plot=TRUE)
dev.off()

png("cluster.png",width=800,height=600)
clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
dev.off()

png("PCA.png",width=800,height=600)
PCASamples(meth,scale=FALSE,center=FALSE)
dev.off()

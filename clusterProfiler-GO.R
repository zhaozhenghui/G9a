library(clusterProfiler)
degenes <- read.csv("GC-CK-up.txt",header = T,stringsAsFactors = F,sep = '\t')
head(degenes)
genelist <- degenes$Gene_symbol
head(genelist)
genelist[duplicated(genelist)]
tr_genelist <- bitr(genelist, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Mm.eg.db")
go <- enrichGO(genelist, OrgDb = org.Mm.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'SYMBOL')
head(go)
write.table(go,"GC-CK-high-methyl.txt")
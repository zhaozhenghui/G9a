

G9a <- read.csv(file = "C:/.../G9a/Counts.csv",sep = ",",header = TRUE)
nrow(G9a)
library(DESeq2)

condition <- factor(c(rep("G9aWT",3),rep("G9aKO",3)),levels = c("G9aWT","G9aKO"))
sampleNames <- c("C1_count","C2_count","C6_count","K3_count","K4_count","K5_count")
colData <- data.frame(name=c("C1_count","C2_count","C6_count","K3_count","K4_count","K5_count"),condition)
row.names(colData) <- sampleNames
dds <- DESeqDataSetFromMatrix(countData=countMatrix,colData=colData,design = ~ condition)
dds <- dds[rowSums(counts(dds))>1,]
nrow(dds)

dds2 <- DESeq(dds)
res <- results(dds2)

resdata <- merge(as.data.frame(res),as.data.frame(counts(dds2,normalize=TRUE)),by="row.names",sort=FALSE)
diff_gene <- subset(resdata, pvalue < 0.05 & abs(log2FoldChange) > 1)
nrow(diff_gene)
up_diff <- subset(resdata,pvalue < 0.05 & log2FoldChange >1)
colnames(up_diff)[1]<-c("gene")
nrow(up_diff)
down_diff <- subset(resdata, log2FoldChange <  -1 & pvalue < 0.05)
colnames(down_diff)[1]<-c("gene")
nrow(down_diff)
write.csv(diff_gene,file = "C:/Users/DELL/Documents/G9a/diff_gene.csv")
write.csv(up_diff,file = "C:/Users/DELL/Documents/G9a/up_diff.csv")
write.csv(down_diff,file = "C:/Users/DELL/Documents/G9a/down_diff.csv")
vsd <- vst(dds2, blind=FALSE)
file.path("vsd", fsep = .Platform$file.sep)



BiocManager::install("DOSE")
library(DOSE)
BiocManager::install("GO.db")
library(GO.db)
BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
BiocManager::install("topGO")
library(topGO)
BiocManager::install("PGSEA")
library(PGSEA)
BiocManager::install("clusterProfiler")
library(clusterProfiler)
up_diff1 <- read.csv("C:/Users/DELL/Documents/G9a/up_diff.csv",header = T)
eg = bitr(up_diff$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
View(eg)
id = as.character(eg[,2])
ego <- enrichGO(gene = id,
                OrgDb = org.Mm.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
View(as.data.frame(ego))
write.csv(ego, file = "C:/Users/DELL/Documents/G9a/up_diffgo.csv")
dotplot(ego)
barplot(ego, showCategory=5)


down_diff1 <- read.csv("C:/Users/DELL/Documents/G9a/down_diff.csv",header = T)
eg = bitr(down_diff1$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
View(eg)
id = as.character(eg[,2])
ego <- enrichGO(gene = id,
                OrgDb = org.Mm.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
View(as.data.frame(ego))
write.csv(ego, file = "C:/Users/DELL/Documents/G9a/down_diffgo.csv")
dotplot(ego)
barplot(ego, showCategory=5)

BiocManager::install("VennDiagram")
library(VennDiagram)
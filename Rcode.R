library(WGCNA)
library(RColorBrewer)
library(GEOquery)
library(org.Hs.eg.db)
data=read.table("GSE256068_processed_data.txt",
                sep = " ",
                fileEncoding = "UTF-16LE",
                fill = T,
                check.names = F)
group=read.table("GSE256068_series_matrix.txt",sep="\t",fill = T)
group=data.frame(t(group))
group=group[-1,]
group=group[,c(19,1,8,13,14)]
group=separate(data = group, col = X13, into = c("X", "type"), sep = ": ")
traits=group[,-4]
View(traits)
data$ENSEMBL=rownames(data)
id=bitr(data$ENSEMBL,toType = "SYMBOL",fromType = "ENSEMBL",OrgDb = org.Hs.eg.db)
datExpr=merge(id,data,by="ENSEMBL")
datExpr = datExpr %>% distinct(SYMBOL, .keep_all = T)
rownames(datExpr)=datExpr[,2]
datExpr=datExpr[,-1:-2]
datExpr=na.omit(datExpr)
m.mad <- apply(datExpr,1,mad)
datExprVar <- datExpr[which(m.mad >max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
datExpr = t(datExprVar[order(apply(datExprVar,1,mad), decreasing = T)[1:8000],])#转置
gsg = goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK##TRUE可用
enableWGCNAThreads()
sampleTree = hclust(dist(datExpr), method = "average")
sizeGrWindow(12,9) #视图大小
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)##聚类检验离群值
sample_colors <- numbers2colors(as.numeric(factor(traits$type)),   colors = rainbow(6),signed = FALSE)
par(mar = c(1,4,3,1),cex=0.8)
plotDendroAndColors(sampleTree, sample_colors,
                    groupLabels = colnames(sample),
                    cex.dendroLabels = 0.8,
                    marAll = c(1, 4, 3, 1),
                    cex.rowText = 0.01,
                    main = "Sample dendrogram and trait heatmap")
##软阈值
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 1.5
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.85,col="red")
##平均连接度
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
###分步
adjacency = adjacency(datExpr, power = 12)
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average");
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
minModuleSize = 30
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs)
METree = hclust(as.dist(MEDiss), method = "average")
MEDissThres = 0.25
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs
# 统计mergedmodule
table(mergedColors)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
ab=data.frame(colnames(datExpr))
ab$color=mergedColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
## module-trait relationship
nSamples = nrow(datExpr)
design=model.matrix(~0+ traits$type)
(colnames(design)=levels(as.factor(traits$type)))
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0); ##不同颜色的模块的ME值矩阵(样本vs模块)
moduleTraitCor = cor(MEs, design , use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
sizeGrWindow(10,6)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar =c(7, 15, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
#####再计算表型与基因的相关性矩阵
recovery = as.data.frame(design[,5]);
names(recovery) = "TLE-HS"
geneTraitSignificance = as.data.frame(cor(datExpr, recovery, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(recovery), sep="");
names(GSPvalue) = paste("p.GS.", names(recovery), sep="");
module = "skyblue"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for TLE-HS",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
####模块与模块之间的相关性
MEs_col = net$MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap",
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T,
                      xLabelsAngle = 90)
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
recovery = as.data.frame(design[,3]);
names(recovery) = "TLE-HS"
MET = orderMEs(cbind(MEs, recovery))
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)
module = "skyblue"
dat=datExpr[,moduleColors==module]
n=t(scale(dat))
group_list=traits$type
ac=data.frame(g=group_list)
rownames(ac)=colnames(n)
colnames(ac) <- "Group"
pheatmap(n, #表达数据
         cluster_rows = T,#行聚类
         cluster_cols = T,#列聚类
         annotation_col =ac,
         annotation_legend=TRUE, # 显示样本分类
         show_rownames = T,# 显示行名
         show_colnames = F,# 显示列名
         scale = "row", #对行标准化
         color =colorRampPalette(c("#0020C2", "#ffffff","#DC143C"))(100),border_color = NA)
dir.create("Module")
setwd("Module")
co=data.frame(table(moduleColors))
data_list <- list()
for (i in co$moduleColors){
  module = i
  probes = colnames(datExpr)
  inModule = (moduleColors==module)
  data_list [[i]] <- data.frame(probes[inModule])
}
out_fileName <- sapply(names(data_list ),function(x){
  paste(x, ".csv", sep='')})
outPath <-getwd()
out_filePath  <- sapply(out_fileName, function(x){
  paste(outPath ,x,sep='/')})
for(i in 1:length(data_list)){
  write.csv(data_list[[i]], file=out_filePath[i])
}
setwd("..")
##批量功能富集
dir.create("Module_KEGG")
setwd("Module_KEGG")
for (i in co$moduleColors){
  module = i
  list_gene=data_list[[i]]
  id=bitr(list_gene[,1],fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
  KEGG<-enrichKEGG(id$ENTREZID,
                   organism = "hsa",
                   pvalueCutoff = 0.05,
                   qvalueCutoff =0.05,
                   minGSSize =30)
  KEGG=setReadable(KEGG,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
  KEGG_list [[i]] <- KEGG@result
}
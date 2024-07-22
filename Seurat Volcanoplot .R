pbmc <- JoinLayers(pbmc )
markers_genes <- FindAllMarkers(pbmc, log2FC.threshold = 0.2, test.use = "wilcox",
                                min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE, max.cells.per.ident = 50,
                                assay = "RNA")
markers_genes %>%
  group_by(cluster) %>%
  top_n(-25, p_val_adj) -> top25

##barplot
mypar(2, 5, mar = c(4, 6, 3, 1))
for (i in unique(top25$cluster)) {
  barplot(sort(setNames(top25$avg_log2FC, top25$gene)[top25$cluster == i], F),
          horiz = T, las = 1, main = paste0(i), border = "white", yaxs = "i")
  abline(v = c(0, 0.25), lty = c(1, 2))
}

##heatmap
markers_genes %>%
  group_by(cluster) %>%
  top_n(-5, p_val_adj) -> top5
pbmc <- ScaleData(pbmc, features = as.character(unique(top5$gene)), assay = "RNA")
DoHeatmap(pbmc,
          features = as.character(unique(top5$gene)), 
          group.by = "cell",
          assay = "RNA")+
  theme(axis.text.y=element_text(size = 13,face="plain",colour = "black"))+
  guides(color=guide_legend(override.aes = list(size=8,alpha=1)))+
  labs(fill = "cell")

##
DotPlot(pbmc, features = rev(as.character(unique(top88$gene))), group.by = "cell",assay = "RNA") + coord_flip()+ylab("Celltype")
###
cell_selection <- subset(pbmc, cells = colnames(pbmc)[pbmc@meta.data[, 'cell'] =="Epithelial"])
cell_selection <- SetIdent(cell_selection, value = "Group")
# Compute differentiall expression
DGE_cell_selection <- FindAllMarkers(cell_selection, log2FC.threshold = 0.2, test.use = "wilcox",
                                     min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE, max.cells.per.ident = 50,
                                     assay = "RNA")
DGE_cell_selection %>%
  group_by(cluster) %>%
  top_n(-5, p_val) -> top5_cell_selection

VlnPlot(cell_selection, features =as.character(unique(top5_cell_selection$gene)),ncol = 3, group.by = "orig.ident", assay = "RNA", pt.size = 0.1)
VlnPlot(cell_selection, features = as.character(unique(top5_cell_selection$gene)),ncol = 5, group.by = "Group", assay = "RNA", pt.size = 0.1,cols = c("#F75D59","#38ACEC"))


##pseudobulking

pseudo_ifnb <- AggregateExpression(scedata, assays = "RNA", return.seurat = T, group.by = c("orig.ident", "new","celltype"))##new是样本ID

pseudo_ifnb$new2=substr(pseudo_ifnb$orig.ident, 1,2)
pseudo_ifnb$celltype.stim <- paste(pseudo_ifnb$celltype, pseudo_ifnb$new2, sep = "_")

Idents(pseudo_ifnb) <- "celltype.stim"

bulk.mono.de <- FindMarkers(object = pseudo_ifnb, 
                            ident.1 = "Oligodendrocyte_IL", 
                            ident.2 = "Oligodendrocyte_CL",
                            test.use = "DESeq2",min.cells.group = 2)
head(bulk.mono.de, n = 15)
bulk.mono.de$gene=rownames(bulk.mono.de)
sig=bulk.mono.de[bulk.mono.de$p_val<0.05,]
sig$cluster=ifelse(sig$avg_log2FC>0,"IL","CL")
sig$gene=rownames(sig)


library(ggplot2)
library(ggrepel)
n=0
bulk.mono.de[which(bulk.mono.de$avg_log2FC >= n & bulk.mono.de$p_val < 0.05),'sig'] <- 'Up'
bulk.mono.de[which(bulk.mono.de$avg_log2FC <= -n & bulk.mono.de$p_val < 0.05),'sig'] <- 'Down'
bulk.mono.de[which(abs(bulk.mono.de$avg_log2FC) <= n | bulk.mono.de$p_val >= 0.05),'sig'] <- 'NoSignifi'
ggplot(bulk.mono.de,aes(x=avg_log2FC,y=-log10(p_val),color=sig))+
  geom_point()+
  scale_color_manual(values=c("#00008B","#808080","#DC143C"))+
  geom_text_repel(
    data = bulk.mono.de[bulk.mono.de$p_val<0.05&abs(bulk.mono.de$avg_log2FC)>n,],
    aes(label = gene),
    size = 5,
    segment.color = "black", show.legend = FALSE )+
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.text.x=element_text(size=14,face="plain",color = "black"), #设置x轴刻度标签的字体属性
        axis.text.y=element_text(size=14,face="plain",color = "black"), #设置y轴刻度标签的字体属性
        axis.title.y=element_text(size = 14,face="plain",color = "black"), #设置y轴的标题的字体属性
        axis.title.x=element_text(size = 14,face="plain",color = "black"))+
  ylab('-log10 (P_Val)')+
  xlab('log2 (FoldChange)')+
  geom_vline(xintercept=c(0),lty=3,col="black",lwd=0.5) +
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)+
  guides(color=guide_legend(override.aes = list(size=8,alpha=1)))



cell_selection_BMYO <- subset(MG1, cells = colnames(MG1)[MG1@meta.data[, "celltype"] == "BMYO"])
cell_selection_BMYO <- SetIdent(cell_selection_BMYO, value = "Group")

DGE_cell_selection<- FindAllMarkers(cell_selection_BMYO, log2FC.threshold = 0.2, test.use = "wilcox", min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE, max.cells.per.ident = 50, assay = "RNA")



DGE_cell_selection %>%
  group_by(cluster) %>%
  top_n(-6, p_val) -> top5_cell_selection
VlnPlot(cell_selection_LHS, features = as.character(unique(top5_cell_selection$gene)),ncol = 6, group.by = "Group", assay = "RNA", pt.size = 0.1,cols = c("#F75D59","#38ACEC"))

write.csv(DGE_cell_selection,"DEG_BMYO.csv")



VlnPlot(cell_selection_LHS, features = as.character(unique(top5_cell_selection$gene)),ncol = 5, group.by = "Group", assay = "RNA", pt.size = 0.1,cols = c("#F75D59","#38ACEC"))




cell_selection_BMYO <- subset(MG1, cells = colnames(MG1)[MG1@meta.data[, "celltype"] == "LASP"])
cell_selection_BMYO <- SetIdent(MG, value = "Group")
DGE_cell_selection3 <- FindAllMarkers(cell_selection_BMYO, log2FC.threshold = -Inf, test.use = "wilcox", min.pct = 0.1, min.diff.pct = 0, only.pos = F, max.cells.per.ident = 50, assay = "RNA")
DGE_cell_selection3=DGE_cell_selection3[DGE_cell_selection3$avg_log2FC>0,]
DGE_cell_selection3$logFC=ifelse(DGE_cell_selection3$cluster=="MT",DGE_cell_selection3$avg_log2FC,-DGE_cell_selection3$avg_log2FC)
rownames(DGE_cell_selection3)=DGE_cell_selection3$gene
write.csv(DGE_cell_selection3,"All gene_Epithelial.csv")

DGE_cell_selection3$diff=abs(DGE_cell_selection3$pct.1 - DGE_cell_selection3$pct.2)

n=0.2
DGE_cell_selection3[which(DGE_cell_selection3$logFC >= n & DGE_cell_selection3$pct.1>0.1&DGE_cell_selection3$p_val<0.01)&DGE_cell_selection3$diff>0.2,'sig'] <- 'Up'
DGE_cell_selection3[which(DGE_cell_selection3$logFC <= -n & DGE_cell_selection3$pct.1>0.1&DGE_cell_selection3$p_val<0.01&DGE_cell_selection3$diff>0.2),'sig'] <- 'Down'
DGE_cell_selection3[which(abs(DGE_cell_selection3$logFC ) <= n | DGE_cell_selection3$pct.1<0.1|DGE_cell_selection3$p_val>0.01|DGE_cell_selection3$diff<0.2),'sig'] <- 'NoSignifi'
diff_gene_deseq2 <- data.frame(subset(DGE_cell_selection3, sig %in% c('Up', 'Down')))


ggplot(DGE_cell_selection3,aes(x=logFC,y=-log10(p_val),color=sig))+
  geom_point(size=5,alpha=0.7)+
  scale_color_manual(values=c("#00008B","#808080","#DC143C"))+
  geom_text_repel(
    data = DGE_cell_selection3[DGE_cell_selection3$p_val<0.01&abs(DGE_cell_selection3$logFC)>n&DGE_cell_selection3$pct.1>0.1&DGE_cell_selection3$diff>0.2,],
    aes(label = gene),
    size = 5,
    segment.color = "black", show.legend = FALSE )+
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.text.x=element_text(size=14,face="plain",color = "black"), #设置x轴刻度标签的字体属性
        axis.text.y=element_text(size=14,face="plain",color = "black"), #设置y轴刻度标签的字体属性
        axis.title.y=element_text(size = 14,face="plain",color = "black"), #设置y轴的标题的字体属性
        axis.title.x=element_text(size = 14,face="plain",color = "black"))+
  ylab('-log10 (P Value)')+
  xlab('log2 (FoldChange)')+
  geom_vline(xintercept=c(-0.2,0.2),lty=3,col="black",lwd=0.5) +
  geom_hline(yintercept = -log10(0.01),lty=3,col="black",lwd=0.5)+
  guides(color=guide_legend(override.aes = list(size=8,alpha=1)))






cluster1.markers <- DGE_cell_selection2 %>%
  mutate(Difference = pct.1 - pct.2) %>%
  rownames_to_column("gene")
log2FC = 1
padj = 0.05
cluster1.markers$threshold="NS";
cluster1.markers[which(cluster1.markers$avg_log2FC  > log2FC & cluster1.markers$p_val <padj),]$threshold="Up";
cluster1.markers[which(cluster1.markers$avg_log2FC  < (-log2FC) & cluster1.markers$p_val < padj),]$threshold="Down";
cluster1.markers$threshold=factor(cluster1.markers$threshold, levels=c('Down','NS','Up'))



ggplot(cluster1.markers, aes(x=Difference, y=avg_log2FC, color = threshold)) +
  geom_point(size=2) +
  scale_color_manual(values=c( "blue","grey","red") ) +
  geom_label_repel(data=subset(cluster1.markers, avg_log2FC >= 1 & Difference >= 0.15 & p_val <= 0.05),aes(label=gene), color="black")+
  geom_label_repel(data=subset(cluster1.markers, avg_log2FC <= -1 & Difference <= -0.15 & p_val <= 0.05),aes(label=gene),color="black")+
  geom_vline(xintercept = 0.0,linetype=2)+
  geom_hline(yintercept = 0,linetype=2)+
  theme_classic()+
  theme(legend.text = element_text(size = 12),#图例大小
        axis.text.x=element_text(colour="black",size=14), #设置x轴刻度标签的字体属性
        axis.text.y=element_text(size=14,face="plain"), #设置y轴刻度标签的字体属性
        axis.title.y=element_text(size = 14,face="plain"), #设置y轴的标题的字体属性
        axis.title.x=element_text(size = 14,face="plain"),
        legend.title =element_text(colour="black",size=14))#设置x轴的标题的字体属性


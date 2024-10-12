library(Seurat)
library(CellChat)
Con=pbmc[,pbmc$Group %in% c( "Control")]
OA=pbmc[,pbmc$Group %in% c( "OA")]
cellChat_Con <- createCellChat(object = Con, group.by = "celltype", assay = "RNA")
cellChat_OA <- createCellChat(object = OA, group.by = "celltype", assay = "RNA")
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search =c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact"), key = c("annotation"))
cellChat_Con@DB <- CellChatDB.use
cellChat_OA@DB <- CellChatDB.use
ptm = Sys.time()
cellChat_Con <- subsetData(cellChat_Con) # This step is necessary even if using the whole database
future::plan("multisession", workers = 8) # do parallel

cellChat_Con <- identifyOverExpressedGenes(cellChat_Con)
cellChat_Con<- identifyOverExpressedInteractions(cellChat_Con)
#> The number of highly variable ligand-receptor pairs used for signaling inference is 692
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))
cellChat_Con <- computeCommunProb(cellChat_Con, type = "triMean")
cellChat_Con <- filterCommunication(cellChat_Con, min.cells = 10)
df.net_Con <- subsetCommunication(cellChat_Con)
cellChat_Con <- computeCommunProbPathway(cellChat_Con)
cellChat_Con <- aggregateNet(cellChat_Con)
ptm = Sys.time()
cellChat_OA <- subsetData(cellChat_OA) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellChat_OA <- identifyOverExpressedGenes(cellChat_OA)
cellChat_OA<- identifyOverExpressedInteractions(cellChat_OA)
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

options(future.globals.maxSize = 1000 * 1024^2)
cellChat_OA <- computeCommunProb(cellChat_OA, type = "triMean")
cellChat_OA <- filterCommunication(cellChat_OA, min.cells = 10)
df.net_OA <- subsetCommunication(cellChat_OA)
cellChat_OA <- computeCommunProbPathway(cellChat_OA)
cellChat_OA <- aggregateNet(cellChat_OA)
cellChat_Con <- netAnalysis_computeCentrality(cellChat_Con, slot.name = "netP")
cellChat_OA <- netAnalysis_computeCentrality(cellChat_OA, slot.name = "netP")
object.list <- list(Con = cellChat_Con, OA = cellChat_OA)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))


pdf("bar.pdf",width = 8,height = 6,pointsize = 21)
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
dev.off()

pdf("netVisual_diffInteraction.pdf",width = 8,height = 6,pointsize = 21)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
gg1 + gg2
dev.off()

pdf("netVisual_diffInteraction.pdf",width = 8,height = 6,pointsize = 21)
gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2
dev.off()

gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

##signalingRole_scatter
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)


##ranNet
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)
gg1 + gg2
##
library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))


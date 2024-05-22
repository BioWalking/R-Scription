library(survival)
library(survminer)
library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)
library(DESeq2)
library(dplyr)
library(TCGAbiolinks)
library(data.table)
library(limma)
library(edgeR)
query <- GDCquery(project = "TCGA-LAML", 
data.category = "Transcriptome Profiling",
data.type = "Gene Expression Quantification",#选定要下载的数据类型
workflow.type = "STAR - Counts"
)
GDCdownload(query, method = "api", files.per.chunk = 100)
expdat <- GDCprepare(query = query)
tpm_matrix=assay(expdat,i="tpm_unstrand")
tpm_matrix=data.frame(tpm_matrix)
tpm_matrix$ENSEMBL=rownames(tpm_matrix)
tpm_matrix$ENSEMBL=unlist(str_split(tpm_matrix$ENSEMBL,"[.]",simplify=T))[,1]
id=bitr(tpm_matrix$ENSEMBL,fromType = "ENSEMBL",toType = c("SYMBOL"),OrgDb = org.Hs.eg.db)
tpm=merge(id,tpm_matrix,by="ENSEMBL",all=FALSE)
tpm<- tpm%>%distinct(`SYMBOL`, .keep_all = T)
rownames(tpm)=tpm[,2]
tpm=tpm[,-1:-2]
rm(tpm_matrix)
clinical <- GDCquery_clinic(project = "TCGA-LAML", type = "clinical")
t_needed=c("submitter_id",
"vital_status",
"days_to_last_follow_up",
"days_to_death")
meta=clinical[,t_needed] #筛选需要的临床信息
meta=meta[meta$vital_status %in% c('Alive','Dead'),]
meta$days_to_last_follow_up[is.na(meta$days_to_last_follow_up)] = 0 
meta$days_to_death[is.na(meta$days_to_death)] = 0
meta$days<-ifelse(meta$vital_status=='Alive',meta$days_to_last_follow_up,meta$days_to_death)
meta$month=round(meta$days/30,2)
ins_genes=rownames(tpm)
design=data.frame(colnames(tpm))
design$submitter_id<- str_sub(design$colnames.tpm.,1,12)
design$submitter_id=gsub('[.]', '-', design$submitter_id)
sur=merge(design,meta,by="submitter_id")
tpm=as.matrix(tpm)
result=data.frame()
for (gene in ins_genes) {
sur$group=ifelse(tpm[gene,]>median(tpm[gene,]),'high','low')
survData = Surv(time=sur$month,
                event=sur$vital_status=='Dead')
KMfit <- survfit(survData ~ sur$group) 
KM_result=data.frame(surv_pvalue(KMfit, method = "survdiff",data = sur))
KM_result[,1]=gene
result=rbind(result,KM_result)
}
sig_genes=result[which(result$pval<0.05),]
colnames(sig_genes)[1]="SYMBOL"
write.csv(sig_genes, file = "sig_genes.csv")
##差异
phe=read.table("GTEX_phenotype.gz",header = T,sep="\t")
RSEM=fread("TcgaTargetGtex_gene_expected_count.gz",check.names = F)
RSEM=data.frame(RSEM)
RSEM$sample=unlist(str_split(RSEM$sample,"[.]",simplify=T))[,1]
colnames(RSEM)[1] <- "ENSEMBL"
id=bitr(RSEM$ENSEMBL,toType = "SYMBOL",fromType = "ENSEMBL",OrgDb = org.Hs.eg.db)
RSEM=merge(id,RSEM,by="ENSEMBL")
RSEM<- RSEM%>%distinct(`SYMBOL`, .keep_all = T)
rownames(RSEM)=RSEM[,2]
RSEM=RSEM[,-1:-2]
sampleid=data.frame(colnames(RSEM))
sampleid <- sampleid[grep("K.562", sampleid$colnames.RSEM.),]
design$sampleid=substr(design$colnames.tpm.,1,15)
pick=c(design$sampleid,sampleid)
pick=data.frame(pick)
sampleid2=data.frame(colnames(RSEM))
colnames(pick) <- "id"
colnames(sampleid2) <- "id"
zz=merge(sampleid2,pick,by="id")
data=RSEM[,zz$id]
list <- c(rep("N", 70), rep("T",147))
DGElist <- DGEList(counts =data, group = list)
DGElist <- calcNormFactors( DGElist )
list <- model.matrix(~factor(list)+0)
colnames(list) <- c("N", "T")
v <- voom(data, plot = TRUE, normalize = "quantile")
ee=data.frame(v$E)
df.fit <- lmFit(ee, list)
df.matrix <- makeContrasts(T - N, levels = list)
fit <- contrasts.fit(df.fit, df.matrix)
fit <- eBayes(fit)
tempOutput <- topTable(fit,n = Inf, adjust = "fdr")
tempOutput$SYMBOL=rownames(tempOutput)
n=2
tempOutput [which(tempOutput $logFC >= n & tempOutput $adj.P.Val< 0.05),'sig'] <- 'Up'
tempOutput [which(tempOutput $logFC <= -n & tempOutput $adj.P.Val < 0.05),'sig'] <- 'Down'
tempOutput [which(abs(tempOutput $logFC) <= n | tempOutput $adj.P.Val >= 0.05),'sig'] <- 'NoSignifi'
All_diffSig <- data.frame(subset(tempOutput , sig %in% c('Up', 'Down')))
write.csv(tempOutput, file = "Allgene.csv")
write.csv(All_diffSig, file = "All_diffSig.csv")
##交集
KM_DEG=merge(sig_genes,All_diffSig,by="SYMBOL")
write.csv(KM_DEG, file = "KM_DEG.csv")



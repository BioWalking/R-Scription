Module_KEGG=data.frame()
for (i in co$moduleColors){
  KEGG=KEGG_list[[i]]
  KEGG$Bg1=as.numeric(unlist(str_split(KEGG$BgRatio,"[/]",simplify=T))[,1])
  KEGG$Bg2=as.numeric(unlist(str_split(KEGG$BgRatio,"[/]",simplify=T))[,2])
  KEGG$Bg=KEGG$Bg1/KEGG$Bg2
  KEGG$Gene2=as.numeric(unlist(str_split(KEGG$GeneRatio,"[/]",simplify=T))[,2])
  KEGG$Gene=KEGG$Count/KEGG$Gene2
  KEGG$Fold_Enrichment=KEGG$Gene/KEGG$Bg
  ww=KEGG[order(KEGG$Fold_Enrichment,decreasing =T),]
  ww=ww[1:5,]
  ww$Module=i
  Module_KEGG=data.frame(rbind(ww,Module_KEGG))
}
fs = list.files('.')
samples <- substr(fs,1,nchar(fs)-9)
Module_KEGG=data.frame()
for (i in fs){
  KEGG=read.csv(i)
  KEGG=KEGG[KEGG$pvalue<0.05,]
  KEGG$Bg1=as.numeric(unlist(str_split(KEGG$BgRatio,"[/]",simplify=T))[,1])
  KEGG$Bg2=as.numeric(unlist(str_split(KEGG$BgRatio,"[/]",simplify=T))[,2])
  KEGG$Bg=KEGG$Bg1/KEGG$Bg2
  KEGG$Gene2=as.numeric(unlist(str_split(KEGG$GeneRatio,"[/]",simplify=T))[,2])
  KEGG$Gene=KEGG$Count/KEGG$Gene2
  KEGG$Fold_Enrichment=KEGG$Gene/KEGG$Bg
  ww=KEGG[order(KEGG$Fold_Enrichment,decreasing =T),]
  ww=ww[1:5,]
  i2<- substr(i,1,nchar(i)-9)
  ww$Module=i2
  Module_KEGG=data.frame(rbind(ww,Module_KEGG))
}
Module_KEGG= Module_KEGG%>%distinct(Description, .keep_all = T)
p <- ggplot(Module_KEGG, aes(x=Fold_Enrichment, fct_reorder(factor(Description), Module)))+
  geom_point(aes(size=Count,color=-1*log10(pvalue)))+
  scale_colour_gradient(low="blue", high="red")+
  labs(
    color=expression(-log10(pvalue)),
    size="Count",
    x="Fold Enrichment")+
  theme_bw()+
  theme(axis.text.x = element_text(size=14,face = "plain",colour = "black"),
        axis.text.y = element_text(size=14,face = "plain",colour = "black"),
        axis.title.x = element_text(size=14,colour = "black"),
        axis.title.y = element_blank(),
        legend.text =  element_text(size = 14,colour = "black"),
        legend.title = element_text(size = 14,colour = "black")
  )+scale_size_continuous(range=c(4,10))+scale_y_discrete(position = "right")
col=cluster$Module
col=as.character(col)
cluster <- Module_KEGG[,c(5,19)] %>% as.data.frame()%>%
  mutate(p="")
p2<- 
  ggplot(cluster,aes(p,Description,fill= col))+
  geom_tile() + 
  scale_y_discrete(position="right") +
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(axis.text.y = element_blank(),
        axis.text.x =element_text(
          angle =90,hjust =0.5,vjust = 0.5))+
  labs(fill = "Module")+scale_fill_manual(values=samples)
p%>%
  insert_left(p2, width = .05)
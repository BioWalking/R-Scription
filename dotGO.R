dotGO <- function(x=GO,showCategory=5,pvalue="p",cutoff="p"){
  GO@result$Bg1=as.numeric(unlist(str_split(GO@result$BgRatio,"[/]",simplify=T))[,1])
  GO@result$Bg2=as.numeric(unlist(str_split(GO@result$BgRatio,"[/]",simplify=T))[,2])
  GO@result$Bg=GO@result$Bg1/GO@result$Bg2
  GO@result$Gene2=as.numeric(unlist(str_split(GO@result$GeneRatio,"[/]",simplify=T))[,2])
  GO@result$Gene=GO@result$Count/GO@result$Gene2
  GO@result$Fold_Enrichment=GO@result$Gene/GO@result$Bg
  GOww<- arrange(GO@result,desc(abs(Fold_Enrichment)))%>%
    group_by(ONTOLOGY)%>%dplyr::slice(1:paste0(showCategory))
  
  if(cutoff=="p"){
    GOww=GOww[GOww$pvalue<0.05,]
  }else if (cutoff=="p.adj"){
    GOww=GOww[GOww$p.adjust<0.05,]
  }else{
    GOww=GOww
  }
  if(pvalue=="p"){
    ggplot(GOww, aes(x=Fold_Enrichment, fct_reorder(factor(Description), Fold_Enrichment)))+
      geom_point(aes(size=Count,color=-1*log10(pvalue)))+
      scale_colour_gradient(low="blue", high="red")+
      labs(
        color=expression(-log10(P)),
        size="Count",
        x="Fold Enrichment")+
      theme_bw()+scale_size_continuous(range=c(4,10))+
      theme(axis.text.x = element_text(size=14,face = "plain",colour = "black"),
            axis.text.y = element_text(size=14,face = "plain",colour = "black"),
            axis.title.x = element_text(size=14,colour = "black"),
            axis.title.y = element_blank(),
            legend.text =  element_text(size = 14,colour = "black"),
            legend.title = element_text(size = 14,colour = "black")
      )+facet_grid(ONTOLOGY ~ .,scale="free")
  }else if (pvalue=="p.adj"){
    ggplot(GOww, aes(x=Fold_Enrichment, fct_reorder(factor(Description), Fold_Enrichment)))+
      geom_point(aes(size=Count,color=-1*log10(p.adjust)))+
      scale_colour_gradient(low="blue", high="red")+
      labs(
        color=expression(-log10(Padj)),
        size="Count",
        x="Fold Enrichment")+
      theme_bw()+scale_size_continuous(range=c(4,10))+
      theme(axis.text.x = element_text(size=14,face = "plain",colour = "black"),
            axis.text.y = element_text(size=14,face = "plain",colour = "black"),
            axis.title.x = element_text(size=14,colour = "black"),
            axis.title.y = element_blank(),
            legend.text =  element_text(size = 14,colour = "black"),
            legend.title = element_text(size = 14,colour = "black")
      )+facet_grid(ONTOLOGY ~ .,scale="free")
    
  }else {
    print("Please check your input! Ensure that the input is 'p',or 'p.adj'  ")
  }
}


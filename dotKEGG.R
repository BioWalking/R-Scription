dotKEGG <- function(x,showCategory=10,pvalue="p",cutoff="p"){
  x@result$Bg1=as.numeric(unlist(str_split(x@result$BgRatio,"[/]",simplify=T))[,1])
  x@result$Bg2=as.numeric(unlist(str_split(x@result$BgRatio,"[/]",simplify=T))[,2])
  x@result$Bg=x@result$Bg1/x@result$Bg2
  x@result$Gene2=as.numeric(unlist(str_split(x@result$GeneRatio,"[/]",simplify=T))[,2])
  x@result$Gene=x@result$Count/x@result$Gene2
  x@result$Fold_Enrichment=x@result$Gene/x@result$Bg
  ww=x@result[order(x@result$Fold_Enrichment,decreasing =T),]
  if(cutoff=="p"){
   ww=ww[ww$pvalue<0.05,]
  }else if (cutoff=="p.adj"){
  ww=ww[ww$p.adjust<0.05,]
  }else{
    ww=ww
  }
  
  if(pvalue=="p"){
  ggplot(ww[1:paste0(showCategory),], aes(x=Fold_Enrichment, fct_reorder(factor(Description), Fold_Enrichment)))+
    geom_point(aes(size=Count,color=-1*log10(pvalue)))+
    scale_colour_gradient(low="blue", high="red")+
    labs(
      color=expression(-log10(P)),
      size="Count",
      x="Fold Enrichment")+
    theme_bw()+
    theme(axis.text.x = element_text(size=14,face = "plain",colour = "black"),
          axis.text.y = element_text(size=14,face = "plain",colour = "black"),
          axis.title.x = element_text(size=14,colour = "black"),
          axis.title.y = element_blank(),
          legend.text =  element_text(size = 14,colour = "black"),
          legend.title = element_text(size = 14,colour = "black")
    )+scale_size_continuous(range=c(4,10))
  }else if (pvalue=="p.adj"){
    ggplot(ww[1:paste0(showCategory),], aes(x=Fold_Enrichment, fct_reorder(factor(Description), Fold_Enrichment)))+
      geom_point(aes(size=Count,color=-1*log10(p.adjust)))+
      scale_colour_gradient(low="blue", high="red")+
      labs(
        color=expression(-log10(P.adj)),
        size="Count",
        x="Fold Enrichment")+
      theme_bw()+
      theme(axis.text.x = element_text(size=14,face = "plain",colour = "black"),
            axis.text.y = element_text(size=14,face = "plain",colour = "black"),
            axis.title.x = element_text(size=14,colour = "black"),
            axis.title.y = element_blank(),
            legend.text =  element_text(size = 14,colour = "black"),
            legend.title = element_text(size = 14,colour = "black")
      )+scale_size_continuous(range=c(4,10))
    
  } else {
   print("Please check your input! Ensure that the input is 'p',or 'p.adj'  ")
 }
}
  

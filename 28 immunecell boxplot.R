
TN=group[group$`her2 receptor status:ch1`=="P",]
re1=re[TN$sample,]
re1$group <-TN$ch
re1$sample <- row.names(re1)
TME_New = melt(re1)
colnames(TME_New)=c("Group","Sample","Celltype","Composition")
if(T){
  mytheme <- theme(plot.title = element_text(size = 10,color="black",hjust = 0.5),
                   axis.title = element_text(size = 10,color ="black"),
                   axis.text = element_text(size= 10,color = "black"),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   axis.text.x = element_text(size=13,angle = 45, hjust = 1 ),
                   axis.title.y =  element_text(size=14,face = "bold"),
                   axis.text.y = element_text(size=13),
                   panel.grid=element_blank(),
                   legend.position = "top",
                   legend.text = element_text(size= 12),
                   legend.title= element_text(size= 12)
  ) }
ggplot(TME_New, aes(x = Celltype, y = Composition))+
  labs(y="Immune infiltration",x= NULL,title = NULL)+
  geom_boxplot(aes(fill = Group),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+
  scale_fill_manual(values = c( "#3BB9FF","#FF7F00" ))+
  theme_classic() + mytheme +
  stat_compare_means(aes(group =  Group),
                     label = "p.signif",
                     method = "wilcox.test",
                     hide.ns = F)
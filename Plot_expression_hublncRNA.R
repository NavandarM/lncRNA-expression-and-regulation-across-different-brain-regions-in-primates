
## Barplot for the expression of hub lncRNAs

library("ggplot2")

Cluster_stats <- read.delim("Statistics_of_Cluster8.txt")

Cluster_stats1 <- Cluster_stats[-13,]
Cluster_stats1$X <- factor(Cluster_stats1$X, levels = unique(levels(Cluster_stats1$X)) )
Cluster_stats1$ClusterName <- factor(  Cluster_stats1$ClusterName, levels = c("C8_to_PCG_B","C8_to_PCG_C","C8_to_PCG_D")   )
Cluster_stats1$X <- factor( Cluster_stats1$X, levels = unique(Cluster_stats1$X)  )

p <- ggplot( Cluster_stats1, aes( x=X, y=Target_Number, fill=Organism )  ) + scale_fill_manual(values = c("#CD6666","#9999CB")) + geom_bar( stat="identity", position=position_dodge(width=0.6), width = 0.5 ) +  facet_grid( rows = vars(ClusterName) ) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p +  theme_bw() +
theme(axis.line = element_line(color='black'), plot.background = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank())


save.image("Statistics_of_Cluster8_2_whiteGrids.RData")
savehistory("Statistics_of_Cluster8_2_whiteGrids.Rhistory")
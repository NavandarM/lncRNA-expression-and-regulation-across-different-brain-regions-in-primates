## Script for heatmap and boxplot visualization of the expression of differentially expressed orthologs lncRNAs from both and either of the species.

## Differentially expressed lncRNAs from Human and Chimpanezee
DiffB <- read.delim("Ortho_diff_both.txt")

# Differentially expressed lncRNAs in Chimpanzee, its Human orthologs showed non-differential expression.
DiffC <- read.delim("Diff_in_Chimp_NonDiff_Human.txt")

# Differentially expressed lncRNAs in Human, its Chimpanzee orthologs showed non-differential expression.
DiffH <- read.delim("Diff_in_human_Nondiff_Chimp.txt")

# Here we merged all the three dataframes.
Combine <- rbind(DiffB, DiffC, DiffH)

# Using 'pheatmap' package to plot heatmap and kmeans clustering for the expression profiles.
library("pheatmap")
Heatmap_input <- Combine[,c(9, 11:26)]
row.names(Heatmap_input) <- Heatmap_input$Human1
Heatmap_input$Human1 <- NULL
Kmeans3 <- pheatmap(Heatmap_input[,c(1:8)], filename = "Heatmap_combine3.pdf", cluster_cols = F, scale = "row", labels_row = "", kmeans_k = 10, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdGy")))(100) )
Clusters <- as.data.frame(Kmeans3$kmeans[[1]])
Heatmap_input$Clusters <- Clusters[ match( row.names(Heatmap_input), row.names(Clusters)  )  ,1]

## Extract the clusters and manually re-arrange the clusters for visualization.

C4 <-  Heatmap_input[ grep("\\b1\\b", Heatmap_input$Clusters), ]
C9 <-  Heatmap_input[ grep("\\b2\\b", Heatmap_input$Clusters), ]
C1 <-  Heatmap_input[ grep("\\b3\\b", Heatmap_input$Clusters), ]
C2 <-  Heatmap_input[ grep("\\b4\\b", Heatmap_input$Clusters), ]
C3 <-  Heatmap_input[ grep("\\b5\\b", Heatmap_input$Clusters), ]
C5 <-  Heatmap_input[ grep("\\b6\\b", Heatmap_input$Clusters), ]
C10 <- Heatmap_input[ grep("\\b7\\b", Heatmap_input$Clusters), ]
C7 <-  Heatmap_input[ grep("\\b8\\b", Heatmap_input$Clusters), ]
C6 <-  Heatmap_input[ grep("\\b9\\b", Heatmap_input$Clusters), ]
C8 <-  Heatmap_input[ grep("\\b10\\b", Heatmap_input$Clusters), ]

Heatmap_input_sort <- rbind( C1, C2, C3, C4, C5, C6, C7, C8, C9, C10 )

C1$Recluster <- NULL
C2$Recluster <- NULL
C3$Recluster <- NULL
C4$Recluster <- NULL
C5$Recluster <- NULL
C6$Recluster <- NULL
C7$Recluster <- NULL
C8$Recluster <- NULL
C9$Recluster <- NULL
C10$Recluster <- NULL
Heatmap_input_sort$Recluster <- NULL

C1$ReCluster <- "C1"
C2$ReCluster <- "C2"
C3$ReCluster <- "C3"
C4$ReCluster <- "C4"
C5$ReCluster <- "C5"
C6$ReCluster <- "C6"
C7$ReCluster <- "C7"
C8$ReCluster <- "C8"
C9$ReCluster <- "C9"
C10$ReCluster <- "C10"
Heatmap_input_sort <- rbind( C1, C2, C3, C4, C5, C6, C7, C8, C9, C10 )

## Final heatmap plotting with 'pheatmap' package

# Heatmap for human
pheatmap(Heatmap_input_sort[,c(1:8)], filename = "Heatmap_HumanbasedClustering_hsa.pdf", cluster_cols = F, scale = "row", labels_row = "", cluster_rows = F , color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdGy")))(100) )

# Heatmap for chimpanzee
pheatmap(Heatmap_input_sort[,c(9:16)], filename = "Heatmap_HumanbasedClustering_pan.pdf", cluster_cols = F, scale = "row", labels_row = "", cluster_rows = F , color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdGy")))(100) )

### Boxplots

# Prepare the dataframe for the boxplot

Tmp1 <- data.frame("IDs"=character(),"Expression"=numeric(), "Region" = character(), Organism= character(),  "Cluster"=integer())
for (i in 1:length(Heatmap_input_sort)){
if ( i <= 8){
Area <- colnames(Heatmap_input_sort)[i]
Tmp<- data.frame("IDs" = row.names(Heatmap_input_sort), "Expression"= Heatmap_input_sort[,i], "Region" = Area, Organism="Human", "Cluster"= Heatmap_input_sort[,18])
Tmp1 = rbind(Tmp1, Tmp)
}
if( i >= 9 & i < 17) {
Area <- colnames(Heatmap_input_sort)[i]
Tmp<- data.frame("IDs" = row.names(Heatmap_input_sort), "Expression"= Heatmap_input_sort[,i], "Region" = Area, Organism="Chimp", "Cluster"= Heatmap_input_sort[,18])
Tmp1 = rbind(Tmp1, Tmp)
}
}
head(Tmp1)
Arrangement <- factor(Tmp1$Region, levels = unique(Tmp1$Region))

# Plotting
p1 <- ggplot(Tmp1, aes(x=Arrangement, y=log2(Expression + 1), fill=Organism, color=Organism)) +
geom_boxplot(outlier.shape = NA, lwd=0.25) + ylim(0,11) +
facet_wrap(~Cluster, nrow=2) +
stat_summary(fun=mean, geom="line", aes( group=Organism), position=position_dodge(1), size= 0.8, color="black" ) +
scale_color_manual(values = c("#9899CB","#CD6666")) + scale_fill_manual(values = c("#9899CB","#CD6666")) + theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p1

save.image("Combined_analysis.RData")
savehistory("Combined_analysis.Rhistory")

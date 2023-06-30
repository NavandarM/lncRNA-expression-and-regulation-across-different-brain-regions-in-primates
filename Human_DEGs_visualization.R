## listed below the differentially expressed genes from any of the brain regions.
## Extracted the information of variance stabilized read counts, normalized read counts of respective differentially expressed genes.
## Heatmap and boxplot visualization of the expression of differentially expressed orthologs lncRNAs from both and either of the species.

DEGs <- read.delim("DEGs_in_any_of_the_brainregion.txt")
vst <- read.delim("VarianceStablizedReadCount_ALL.txt")
nrc <- read.delim("NormalizedReadCount_All_human.txt")
DEGs_nrc <- nrc[match(DEGs$GeneIds, nrc$Gene), ]
DEGs_vst <- vst[match(DEGs$GeneIds, row.names(vst)), ]

# Export the variance stabilized read counts, normalized read counts for differentially expressed genes.
write.table(DEGs_nrc,"DEGs_nrc_Human.txt", sep="\t", quote = F)
write.table(DEGs_vst,"DEGs_vst_Human.txt", sep="\t", quote = F)

# Read the normalized read counts for differentially expressed genes
DEGs_nrc <-read.delim("DEGs_nrc_Human.txt")
row.names(DEGs_nrc) <- DEGs_nrc$Gene
DEGs_nrc$Gene <- NULL

# Row Mean the counts of the same brain region.
DEGs_nrc$CB <- rowMeans(DEGs_nrc[,c(6:10)])
DEGs_nrc$STR <- rowMeans(DEGs_nrc[,c(23:27)])
DEGs_nrc$HIP <- rowMeans(DEGs_nrc[,c(33:36)])
DEGs_nrc$ACC <- rowMeans(DEGs_nrc[,c(1:5)])
DEGs_nrc$DPFC <- rowMeans(DEGs_nrc[,c(11:16)])
DEGs_nrc$VPFC <- rowMeans(DEGs_nrc[,c(37:40)])
DEGs_nrc$PMC <- rowMeans(DEGs_nrc[,c(17:22)])
DEGs_nrc$V1C <- rowMeans(DEGs_nrc[,c(28:32)])
DEGs_nrc1 <- DEGs_nrc[,c(41:48)]

# Plot the initial heatmap using pheatmap package and perform k means clustering.
library("pheatmap")
X <- pheatmap(DEGs_nrc1, filename = "DEGs_pcg.kmeans.pdf", kmeans_k = 8, cluster_cols = F, scale = "row")
Cluster <- as.data.frame(X$kmeans[1])
DEGs_nrc$Clusters <- Cluster[ match( row.names(DEGs_nrc), row.names(Cluster) ),1]

c1 <- DEGs_nrc[ grep("\\b1\\b", DEGs_nrc$Clusters), ]
c2 <- DEGs_nrc[ grep("\\b2\\b", DEGs_nrc$Clusters), ]
c3 <- DEGs_nrc[ grep("\\b3\\b", DEGs_nrc$Clusters), ]
c4 <- DEGs_nrc[ grep("\\b4\\b", DEGs_nrc$Clusters), ]
c5 <- DEGs_nrc[ grep("\\b5\\b", DEGs_nrc$Clusters), ]
c6 <- DEGs_nrc[ grep("\\b6\\b", DEGs_nrc$Clusters), ]
c7 <- DEGs_nrc[ grep("\\b7\\b", DEGs_nrc$Clusters), ]
c8 <- DEGs_nrc[ grep("\\b8\\b", DEGs_nrc$Clusters), ]

# Re-arragned the clusters mannualy as per its expression in different brain regions.
c1$Recluster <- "C1"
c5$Recluster <- "C2"
c3$Recluster <- "C3"
c4$Recluster <- "C4"
c2$Recluster <- "C5"
c6$Recluster <- "C6"
c7$Recluster <- "C8"
c8$Recluster <- "C7"
Combined_clusters <- rbind(c1,c5,c3,c4,c2,c6,c8,c7)


# Use the 'dplyr' and 'pheatmap' packages together to get refined version of the heatmap.
library("dplyr")
annotateDF <- select(Combined_clusters, Recluster)
pheatmap(Combined_clusters[,c(41:48)], filename = "DEGs_hsa_pcg_new.pdf", cluster_cols = F, scale = "row", cluster_rows = F, show_rownames = F, annotation_row = annotateDF)
write.table(Combined_clusters, "DEGs_hsa.txt", sep="\t", quote = F)
write.table(Combined_clusters[,c(51,50)], "HSA_order.txt", sep="\t", quote = F)


### Box Plot (Supplementary Figure)

Input_boxPlot <- Combined_clusters[,c(41:48,50)]
head(Input_boxPlot)
H_C <- Input_boxPlot
head(H_C)
Tmp1 <- data.frame("IDs"=character(),"Expression"=numeric(), "Region" = character(), Organism= character(),  "Cluster"=character())
for (i in 1:length(H_C)){
if ( i <= 8){
Area <- colnames(H_C)[i]
Tmp<- data.frame("IDs" = row.names(H_C), "Expression"= H_C[,i], "Region" = Area, Organism="Human", "Cluster"= H_C[,9])
Tmp1 = rbind(Tmp1, Tmp)
}
}
Arrangement <- factor(Tmp1$Region, levels = unique(Tmp1$Region))
Arrangement_clusters <- factor(Tmp1$Cluster, levels = c("C1","C2","C3","C4","C5","C6","C7","C8"))

p1 <- ggplot(Tmp1, aes(x=Arrangement, y=log2(Expression + 1), fill=Organism)) +
geom_boxplot(outlier.shape = NA, lwd=0.25) + ylim(0,16.5) +
facet_wrap(~Cluster, nrow=2) +
stat_summary(fun=mean, geom="line", aes(group=Organism), position=position_dodge(1), size= 0.8 ) +
scale_fill_manual(values = c("#CD6666")) + theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p1








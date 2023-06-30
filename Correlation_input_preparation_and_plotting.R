
library(circlize)
library(chorddiag)  #devtools::install_github("mattflor/chorddiag")

circos.clear()

c1_lnc_mrna <- na.omit(read.delim("C1_lnc_mrna.txt"))
c2_lnc_mrna <- na.omit(read.delim("C2_lnc_mrna.txt"))
c3_lnc_mrna <- na.omit(read.delim("C3_lnc_mrna.txt"))
c4_lnc_mrna <- na.omit(read.delim("C4_lnc_mrna.txt"))
c5_lnc_mrna <- na.omit(read.delim("C5_lnc_mrna.txt"))
c6_lnc_mrna <- na.omit(read.delim("C6_lnc_mrna.txt"))
c7_lnc_mrna <- na.omit(read.delim("C7_lnc_mrna.txt"))
c8_lnc_mrna <- na.omit(read.delim("C8_lnc_mrna.txt"))
c9_lnc_mrna <- na.omit(read.delim("C9_lnc_mrna.txt"))
c10_lnc_mrna <- na.omit(read.delim("C10_lnc_mrna.txt"))

up_curated_c1 <- c1_lnc_mrna[ c1_lnc_mrna$V3 > 0.8 & c1_lnc_mrna$V4 < 0.01, ]
up_curated_c2 <- c2_lnc_mrna[ c2_lnc_mrna$V3 > 0.8 & c2_lnc_mrna$V4 < 0.01, ]
up_curated_c3 <- c3_lnc_mrna[ c3_lnc_mrna$V3 > 0.8 & c3_lnc_mrna$V4 < 0.01, ]
up_curated_c4 <- c4_lnc_mrna[ c4_lnc_mrna$V3 > 0.8 & c4_lnc_mrna$V4 < 0.01, ]
up_curated_c5 <- c5_lnc_mrna[ c5_lnc_mrna$V3 > 0.8 & c5_lnc_mrna$V4 < 0.01, ]
up_curated_c6 <- c6_lnc_mrna[ c6_lnc_mrna$V3 > 0.8 & c6_lnc_mrna$V4 < 0.01, ]
up_curated_c7 <- c7_lnc_mrna[ c7_lnc_mrna$V3 > 0.8 & c7_lnc_mrna$V4 < 0.01, ]
up_curated_c8 <- c8_lnc_mrna[ c8_lnc_mrna$V3 > 0.8 & c8_lnc_mrna$V4 < 0.01, ]
up_curated_c9 <- c9_lnc_mrna[ c9_lnc_mrna$V3 > 0.8 & c9_lnc_mrna$V4 < 0.01, ]
up_curated_c10 <- c10_lnc_mrna[ c10_lnc_mrna$V3 > 0.8 & c10_lnc_mrna$V4 < 0.01, ]

Down_curated_c1 <- c1_lnc_mrna[ c1_lnc_mrna$V3 < -0.8 & c1_lnc_mrna$V4 < 0.01, ]
Down_curated_c2 <- c2_lnc_mrna[ c2_lnc_mrna$V3 < -0.8 & c2_lnc_mrna$V4 < 0.01, ]
Down_curated_c3 <- c3_lnc_mrna[ c3_lnc_mrna$V3 < -0.8 & c3_lnc_mrna$V4 < 0.01, ]
Down_curated_c4 <- c4_lnc_mrna[ c4_lnc_mrna$V3 < -0.8 & c4_lnc_mrna$V4 < 0.01, ]
Down_curated_c5 <- c5_lnc_mrna[ c5_lnc_mrna$V3 < -0.8 & c5_lnc_mrna$V4 < 0.01, ]
Down_curated_c6 <- c6_lnc_mrna[ c6_lnc_mrna$V3 < -0.8 & c6_lnc_mrna$V4 < 0.01, ]
Down_curated_c7 <- c7_lnc_mrna[ c7_lnc_mrna$V3 < -0.8 & c7_lnc_mrna$V4 < 0.01, ]
Down_curated_c8 <- c8_lnc_mrna[ c8_lnc_mrna$V3 < -0.8 & c8_lnc_mrna$V4 < 0.01, ]
Down_curated_c9 <- c9_lnc_mrna[ c9_lnc_mrna$V3 < -0.8 & c9_lnc_mrna$V4 < 0.01, ]
Down_curated_c10 <- c10_lnc_mrna[ c10_lnc_mrna$V3 < -0.8 & c10_lnc_mrna$V4 < 0.01, ]

Stats_c1 <- as.data.frame(table(up_curated_c1$V6))
Stats_c1$Var1 <- paste0(Stats_c1$Var1, "_pcg")
Stats_c1$lncClust <- "C1_lnc"
Stats_c1
Stats_c2 <- as.data.frame(table(up_curated_c2$V6))
Stats_c2$Var1 <- paste0(Stats_c2$Var1, "_pcg")
Stats_c2$lncClust <- "C2_lnc"
Stats_c2
Stats_c3 <- as.data.frame(table(up_curated_c3$V6))
Stats_c3$Var1 <- paste0(Stats_c3$Var1, "_pcg")
Stats_c3$lncClust <- "C3_lnc"
Stats_c3
Stats_c4 <- as.data.frame(table(up_curated_c4$V6))
Stats_c4$Var1 <- paste0(Stats_c4$Var1, "_pcg")
Stats_c4$lncClust <- "C4_lnc"
Stats_c4
Stats_c5 <- as.data.frame(table(up_curated_c5$V6))
Stats_c5$Var1 <- paste0(Stats_c5$Var1, "_pcg")
Stats_c5$lncClust <- "C5_lnc"
Stats_c5
Stats_c6 <- as.data.frame(table(up_curated_c6$V6))
Stats_c6$Var1 <- paste0(Stats_c6$Var1, "_pcg")
Stats_c6$lncClust <- "C6_lnc"
Stats_c6
Stats_c7 <- as.data.frame(table(up_curated_c7$V6))
Stats_c7$Var1 <- paste0(Stats_c7$Var1, "_pcg")
Stats_c7$lncClust <- "C7_lnc"
Stats_c7
Stats_c8 <- as.data.frame(table(up_curated_c8$V6))
Stats_c8$Var1 <- paste0(Stats_c8$Var1, "_pcg")
Stats_c8$lncClust <- "C8_lnc"
Stats_c8
Stats_c9 <- as.data.frame(table(up_curated_c9$V6))
Stats_c9$Var1 <- paste0(Stats_c9$Var1, "_pcg")
Stats_c9$lncClust <- "C9_lnc"
Stats_c9
Stats_c10 <- as.data.frame(table(up_curated_c10$V6))
Stats_c10$Var1 <- paste0(Stats_c10$Var1, "_pcg")
Stats_c10$lncClust <- "C10_lnc"
Stats_c10

Reference_Stats <- read.delim("Reference_lnc_pcg_heatmap_counts.txt")
PCor_Stats_cMaster <- rbind(Stats_c1, Stats_c2, Stats_c4, Stats_c5, Stats_c6, Stats_c7, Stats_c8, Stats_c9, Stats_c10 )
PCor_Stat_cMaster1  <- PCor_Stats_cMaster[,c(3,1,2)]
PCor_Stat_cMaster1$Ref_lncs <- Reference_Stats[ match( PCor_Stat_cMaster1$lncClust, Reference_Stats$Cluster_lnc ), "numbers_lnc"  ]
PCor_Stat_cMaster1$Ref_pcg <- Reference_Stats[ match( PCor_Stat_cMaster1$Var1, Reference_Stats$Cluster_pcg ), "numbers_pcg" ]
PCor_Stat_cMaster1$Multiply <- PCor_Stat_cMaster1$Ref_lncs * PCor_Stat_cMaster1$Ref_pcg
PCor_Stat_cMaster1$Fraction <- PCor_Stat_cMaster1$Freq  / PCor_Stat_cMaster1$Multiply

unique(PCor_Stat_cMaster1$lncClust)
PCor_Stat_cMaster1

unique(PCor_Stat_cMaster1$lncClust)

unique(PCor_Stat_cMaster1$Var1)

mycolor <- c("#FF5733","#B2BD0E","#0EBD72","#0EBDB0","#0E72BD","#D9A1CB","#FFF633","#A82104","#1E6482", rep("#9D9D9D",8))

colnames(PCor_Stat_cMaster1)

circos.par(start.degree = 90, gap.degree = 4, track.margin = c(-0.1, 0.1), points.overflow.warning = FALSE)
par(mar = rep(0, 4))
Order_of_ClustName
unique(PCor_Stat_cMaster1$lncClust)
Order_of_ClustName <- c("C1_lnc" ,"C2_lnc", "C4_lnc","C5_lnc", "C6_lnc", "C7_lnc" ,"C8_lnc","C9_lnc","C10_lnc", "C1_pcg", "C2_pcg","C3_pcg","C4_pcg", "C5_pcg", "C6_pcg", "C7_pcg", "C8_pcg")
mycolor <- c("#FF5733","#B2BD0E","#0EBD72","#0EBDB0","#D9A1CB","#0E72BD","#32B423","#A82104","#1E6482", rep("#9D9D9D",8))
pdf("Postive_Correlation_chordDia_fraction_hsa.pdf")
chordDiagram(
x = PCor_Stat_cMaster1[,c(1,2,7)],
order = Order_of_ClustName,
grid.col = mycolor,
transparency = 0.40,
directional = 1,
direction.type = c("arrows", "diffHeight"),
diffHeight  = -0.04,
annotationTrack = "grid",
annotationTrackHeight = c(0.05, 0.1),
link.arr.type = "big.arrow",
link.sort = TRUE,
link.largest.ontop = TRUE)
circos.trackPlotRegion(
track.index = 1,
bg.border = NA,
panel.fun = function(x, y) {
xlim = get.cell.meta.data("xlim")
sector.index = get.cell.meta.data("sector.index")
# Add names to the sector.
circos.text(
x = mean(xlim),
y = 3.2,
labels = sector.index,
facing = "bending",
cex = 0.8
)
}
)
dev.off()

Stats_c1 <- as.data.frame(table(Down_curated_c1$V6))
Stats_c1$Var1 <- paste0(Stats_c1$Var1, "_pcg")
Stats_c1$lncClust <- "C1_lnc"
Stats_c1
Stats_c2 <- as.data.frame(table(Down_curated_c2$V6))
Stats_c2$Var1 <- paste0(Stats_c2$Var1, "_pcg")
Stats_c2$lncClust <- "C2_lnc"
Stats_c2
Stats_c3 <- as.data.frame(table(Down_curated_c3$V6))
Stats_c3$Var1 <- paste0(Stats_c3$Var1, "_pcg")
Stats_c3$lncClust <- "C3_lnc"
Stats_c3
Stats_c4 <- as.data.frame(table(Down_curated_c4$V6))
Stats_c4$Var1 <- paste0(Stats_c4$Var1, "_pcg")
Stats_c4$lncClust <- "C4_lnc"
Stats_c4
Stats_c5 <- as.data.frame(table(Down_curated_c5$V6))
Stats_c5$Var1 <- paste0(Stats_c5$Var1, "_pcg")
Stats_c5$lncClust <- "C5_lnc"
Stats_c5
Stats_c6 <- as.data.frame(table(Down_curated_c6$V6))
Stats_c6$Var1 <- paste0(Stats_c6$Var1, "_pcg")
Stats_c6$lncClust <- "C6_lnc"
Stats_c6
Stats_c7 <- as.data.frame(table(Down_curated_c7$V6))
Stats_c7$Var1 <- paste0(Stats_c7$Var1, "_pcg")
Stats_c7$lncClust <- "C7_lnc"
Stats_c7
Stats_c8 <- as.data.frame(table(Down_curated_c8$V6))
Stats_c8$Var1 <- paste0(Stats_c8$Var1, "_pcg")
Stats_c8$lncClust <- "C8_lnc"
Stats_c8
Stats_c9 <- as.data.frame(table(Down_curated_c9$V6))
Stats_c9$Var1 <- paste0(Stats_c9$Var1, "_pcg")
Stats_c9$lncClust <- "C9_lnc"
Stats_c9
Stats_c10 <- as.data.frame(table(Down_curated_c10$V6))
Stats_c10$Var1 <- paste0(Stats_c10$Var1, "_pcg")
Stats_c10$lncClust <- "C10_lnc"
Stats_c10
NCor_Stat_cMaster <- rbind(Stats_c1,Stats_c4,Stats_c5,Stats_c6,Stats_c7,Stats_c8, Stats_c9, Stats_c10 )
unique(NCor_Stat_cMaster$lncClust)
unique(NCor_Stat_cMaster$Var1)
NCor_Stat_cMaster1  <- NCor_Stat_cMaster[,c(3,1,2)]
NCor_Stat_cMaster1
NCor_Stat_cMaster1$Ref_lncs <- Reference_Stats[ match( NCor_Stat_cMaster1$lncClust, Reference_Stats$Cluster_lnc ), "numbers_lnc"    ]
NCor_Stat_cMaster1$Ref_pcg <- Reference_Stats[ match( NCor_Stat_cMaster1$Var1, Reference_Stats$Cluster_pcg ), "numbers_pcg"    ]
NCor_Stat_cMaster1
NCor_Stat_cMaster1$Multiply <- NCor_Stat_cMaster1$Ref_lncs * NCor_Stat_cMaster1$Ref_pcg
NCor_Stat_cMaster1$Fraction <-  NCor_Stat_cMaster1$Freq /NCor_Stat_cMaster1$Multiply
NCor_Stat_cMaster1
unique(NCor_Stat_cMaster$lncClust)
unique(NCor_Stat_cMaster$Var1)
Order_of_ClustName <- c("C1_lnc" , "C4_lnc","C5_lnc", "C6_lnc", "C7_lnc" ,"C8_lnc","C9_lnc","C10_lnc", "C1_pcg", "C2_pcg","C3_pcg","C4_pcg", "C5_pcg", "C6_pcg", "C7_pcg", "C8_pcg")
mycolor <- c("#FF5733","#0EBD72","#0EBDB0","#D9A1CB","#0E72BD","#32B423","#A82104","#1E6482", rep("#9D9D9D",8))

pdf("Negative_Correlation_chordDia_fraction_hsa.pdf")
chordDiagram(
x = NCor_Stat_cMaster1[,c(1,2,7)],
order = Order_of_ClustName,
grid.col = mycolor,
transparency = 0.40,
directional = 1,
direction.type = c("arrows", "diffHeight"),
diffHeight  = -0.04,
annotationTrack = "grid",
annotationTrackHeight = c(0.05, 0.1),
link.arr.type = "big.arrow",
link.sort = TRUE,
link.largest.ontop = TRUE)
circos.trackPlotRegion(
track.index = 1,
bg.border = NA,
panel.fun = function(x, y) {
xlim = get.cell.meta.data("xlim")
sector.index = get.cell.meta.data("sector.index")
# Add names to the sector.
circos.text(
x = mean(xlim),
y = 3.2,
labels = sector.index,
facing = "bending",
cex = 0.8
)
}
)
dev.off()



write.table(up_curated_c1, "up_c1_hsa.txt", sep="\t", quote = F, row.names = F)
write.table(up_curated_c2, "up_c2_hsa.txt", sep="\t", quote = F, row.names = F)
write.table(up_curated_c3, "up_c3_hsa.txt", sep="\t", quote = F, row.names = F)
write.table(up_curated_c4, "up_c4_hsa.txt", sep="\t", quote = F, row.names = F)
write.table(up_curated_c5, "up_c5_hsa.txt", sep="\t", quote = F, row.names = F)
write.table(up_curated_c6, "up_c6_hsa.txt", sep="\t", quote = F, row.names = F)
write.table(up_curated_c7, "up_c7_hsa.txt", sep="\t", quote = F, row.names = F)
write.table(up_curated_c8, "up_c8_hsa.txt", sep="\t", quote = F, row.names = F)
write.table(up_curated_c9, "up_c9_hsa.txt", sep="\t", quote = F, row.names = F)
write.table(up_curated_c10, "up_c10_hsa.txt", sep="\t", quote = F, row.names = F)

write.table(Down_curated_c1, "Down_c1_hsa.txt", sep="\t", quote = F, row.names = F)
write.table(Down_curated_c2, "Down_c2_hsa.txt", sep="\t", quote = F, row.names = F)
write.table(Down_curated_c3, "Down_c3_hsa.txt", sep="\t", quote = F, row.names = F)
write.table(Down_curated_c4, "Down_c4_hsa.txt", sep="\t", quote = F, row.names = F)
write.table(Down_curated_c5, "Down_c5_hsa.txt", sep="\t", quote = F, row.names = F)
write.table(Down_curated_c6, "Down_c6_hsa.txt", sep="\t", quote = F, row.names = F)
write.table(Down_curated_c7, "Down_c7_hsa.txt", sep="\t", quote = F, row.names = F)
write.table(Down_curated_c8, "Down_c8_hsa.txt", sep="\t", quote = F, row.names = F)
write.table(Down_curated_c9, "Down_c9_hsa.txt", sep="\t", quote = F, row.names = F)
write.table(Down_curated_c10, "Down_c10_hsa.txt", sep="\t", quote = F, row.names = F)

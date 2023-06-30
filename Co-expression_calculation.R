### Calculate correlation of expression between Human lncRNAs from each clusters and protein coding genes.

File_list = c("Uniq_C1_hsa.txt","Uniq_C2_hsa.txt","Uniq_C3_hsa.txt","Uniq_C4_hsa.txt","Uniq_C5_hsa.txt","Uniq_C6_hsa.txt","Uniq_C7_hsa.txt","Uniq_C8_hsa.txt","Uniq_C9_hsa.txt","Uniq_C10_hsa.txt")
PCG <- read.delim("PCG_vst_for_corHSA.txt")
count = 1

for ( entity in File_list) {
 
Lnc <- read.delim(entity)
C_lnc_mrna <- matrix(ncol = 6)

for(i in 1:nrow(Lnc)){
          for(j in 1:nrow(PCG)){
                   Value <- cor.test( as.numeric(Lnc[i,c(1:40)]), as.numeric(PCG[j,c(1:40)]) )
  C_lnc_mrna <- rbind( C_lnc_mrna, matrix(c( row.names(Lnc)[i], row.names(PCG)[j], Value$estimate , Value$p.value,as.character(Lnc$ReCluster[i]),as.character(PCG$ReCluster[j])), ncol = 6) )
  }
}

New_file= paste0("C",count,"_lnc_mrna.txt")
write.table(as.data.frame(C_lnc_mrna),New_file, sep="\t", quote=F)
}

### Calculate correlation of expression between Chimpanzee lncRNAs from each clusters and protein coding genes.

File_list = c("Uniq_C1_pan.txt", "Uniq_C2_pan.txt", "Uniq_C3_pan.txt", "Uniq_C4_pan.txt", "Uniq_C5_pan.txt", "Uniq_C6_pan.txt", "Uniq_C7_pan.txt", "Uniq_C8_pan.txt", "Uniq_C9_pan.txt","Uniq_C10_pan.txt")

PCG <- read.delim("Chimp_PCG.txt")
count = 1

for ( entity in File_list) {
 
Lnc <- read.delim(entity)
C_lnc_mrna <- matrix(ncol = 6)

for(i in 1:nrow(Lnc)){
          for(j in 1:nrow(PCG)){
                   Value <- cor.test( as.numeric(Lnc[i,c(1:40)]), as.numeric(PCG[j,c(1:40)]) )
  C_lnc_mrna <- rbind( C_lnc_mrna, matrix(c( row.names(Lnc)[i], row.names(PCG)[j], Value$estimate , Value$p.value,as.character(Lnc$ReCluster[i]),as.character(PCG$ReCluster[j])), ncol = 6) )
  }
}

New_file= paste0("C",count,"_lnc_mrna.txt")
write.table(as.data.frame(C_lnc_mrna),New_file, sep="\t", quote=F)

}


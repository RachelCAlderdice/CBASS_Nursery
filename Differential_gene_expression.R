library("DESeq2")
library("ggplot2")

############### Overview of gene expression using DESeq2 ###################

#input StringTie gene count matrix & metadata
countData <- as.matrix(read.csv("Inputfiles/gene_count_matrix_ahyc.csv", row.names= "gene_id"))
meta <- read.csv("Inputfiles/Metadata_CBASS_GBR_ahyc.csv", sep=",", row.names=1)
meta$Temp_Year <- paste(meta$Temperature, meta$Year, sep="_")

#remove samples not passing QC (manually remove any single samples)
colData <- subset(meta, meta$Temperature == "30" |meta$Temperature == "34")

#check names match in files
all(rownames(colData) %in% colnames(countData))
all(rownames(colData) == colnames(countData))

dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData, design = ~ Temp_Year)
dds <- DESeq(dds)

############### PCA analysis of counts ###################

vst=varianceStabilizingTransformation(dds)

pca_data<-plotPCA(vst, intgroup= 'Temp_Year', returnData = TRUE)
pca_data$Colony<-colData$Colony
pca_data$Year <- as.factor(colData$Year)

COL2=c("lightblue", "blue","#FFFFCC","yellow")
percentVar <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data,(aes(x=PC1,y=PC2, group= Year)))+ 
  geom_point(aes(shape=Colony, fill=group), size=5)+
  theme_classic()+
  labs(col="Temp_Year")+
  scale_shape_manual(values=c(21,23,22,24,25)) +
  scale_fill_manual(values=COL2)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  theme(axis.text = element_text(size=20))+
  theme(axis.title = element_text(size=20))+
  theme(legend.title = element_text(size=20))+
  theme(text = element_text(size=20))+
  theme(aspect.ratio = 1)

#ggsave("PCA_plot.pdf", dpi=300,width = 6, height = 4)

############### Differential gene expression using DESeq2 ###################

res1 <-results(dds, contrast=c("Temp_Year", "30_21", "30_22"), lfcThreshold = 0.0, alpha=0.05)
message(summary(res1, alpha=0.05))
results1 <- as.data.frame(subset(res1, res1$padj<0.05))
write.table(results1, "./Results/DE_output_ahya_30_21v30_22_FC0.0_p0.05.txt", sep= "\t", quote = FALSE, row.names = TRUE)

res2 <-results(dds, contrast=c("Temp_Year", "34_21", "34_22"), lfcThreshold = 0.0, alpha=0.05)
message(summary(res2, alpha=0.05))
results2 <- as.data.frame(subset(res2, res2$padj<0.05))
write.table(results2, "./Results/DE_output_ahya_34_21v34_22_FC0.0_p0.05.txt", sep= "\t", quote = FALSE, row.names = TRUE)

res3 <-results(dds, contrast=c("Temp_Year", "30_21", "34_21"), lfcThreshold = 0.0, alpha=0.05)
message(summary(res3, alpha=0.05))
results3 <- as.data.frame(subset(res3, res3$padj<0.05))
write.table(results3, "./Results/DE_output_ahya_30_21v34_21_FC0.0_p0.05.txt", sep= "\t", quote = FALSE, row.names = TRUE)

res4 <-results(dds, contrast=c("Temp_Year", "30_22", "34_22"), lfcThreshold = 0.0, alpha=0.05)
message(summary(res4, alpha=0.05))
results4 <- as.data.frame(subset(res4, res4$padj<0.05))
write.table(results4, "./Results/DE_output_ahya_30_22v34_22_FC0.0_p0.05.txt", sep= "\t", quote = FALSE, row.names = TRUE)


################################### Barplot of DEG regulation ############################################

#input DEG list 
DE<-read.table("Results/DE_output_ahya_30_21v30_22_FC0.0_p0.05.txt", row.names = 1, header = TRUE) #change accordingly
DE_up=subset(DE, DE$log2FoldChange > 0) 
DE_down=subset(DE, DE$log2FoldChange < 0) 

#record values & input table
DE<- read.delim("Inputfiles/DEG_regulation_barplot.txt", sep="\t")
DE$Regulation <- as.character(DE$Regulation)
DE$Group <- factor(DE$Group, levels=unique(DE$Group))

ggplot(data = DE, aes(x = Group, y = genes, fill= Regulation))+
  ylim(-1600, 2000)+
  geom_col(colour= "black", width= 0.8, size=0.9)+
  scale_x_discrete(breaks=c("30_21v30_22","34_21v34_22", "30_21v34_21", "30_22v34_22")) +
  theme_classic()+
  scale_fill_manual(values=c("#666666", "#FFFFFF"))+
  theme(aspect.ratio = 1)+
  ylab(label="Genes")+
  xlab(label= "Group")+
  theme(legend.text=element_text(size=18))+
  theme(axis.text = element_text(size=18))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5))+
  theme(axis.title = element_text(size=22))+
  theme(text = element_text(size=20)) +
  theme(axis.line = element_line(colour = 'black', size = 1.5))+
  theme(axis.ticks = element_line(colour = "black", size = 1.5)) +
  theme(axis.ticks.length=unit(.25, "cm"))

#ggsave("Barplot_DEG_regulation.pdf", dpi=300,width = 6, height = 8)
#finalise labels in affinitydesigner
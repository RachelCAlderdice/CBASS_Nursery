library("topGO")
library("data.table")
library("ggplot2")
library("RColorBrewer")
library("ggVennDiagram")
library("dplyr")

############### GO enrichment using topGO ###################

#select DEGs from eggnog annotation file
annot=read.delim("Inputfile/MM_g9x4xoep.emapper.annotations.tsv")
annot$X.query=gsub("\\.t1", "", annot$X.query)
annot$X.query=gsub("\\.t2", "", annot$X.query)
Allgenes=annot[c("X.query", "GOs")]
Allgenes=subset(Allgenes, !(Allgenes$GOs == "-"))
write.table(Allgenes, "Allgenes_GO_annot_new.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

geneID2GO = readMappings(file = "Allgenes_GO_annot_new.txt") #gene names and GO list in 2 columns of all samples
#all possible gene names with GO annotation
geneUniverse = names(geneID2GO)
message("Number of genes with GO annotation: ", length(geneUniverse)) 

DE <- read.delim("Results/DE_output_ahya_34_21v34_22_FC0.0_p0.05.txt", sep="\t", header= TRUE) #change DEG list accordingly

results1$gene_names <- row.names(results1)
DE_genes <- results1$gene_names

keep=DE_genes %in% geneUniverse
keep=which(keep==TRUE)
DE_genes=DE_genes[keep]
#make named vector list of factors showing which are GOI
geneList=factor(as.integer(geneUniverse %in% DE_genes))
names(geneList)=geneUniverse
message("Number of DE genes with GO annotation: ", length(intersect(geneUniverse,DE_genes)))

DE$gene_names <- row.names(DE)
DE_genes <- DE$gene_names

keep=DE_genes %in% geneUniverse
keep=which(keep==TRUE)
DE_genes=DE_genes[keep]
#make named vector list of factors showing which are GOI
geneList=factor(as.integer(geneUniverse %in% DE_genes))
names(geneList)=geneUniverse
message("Number of DE genes with GO annotation: ", length(intersect(geneUniverse,DE_genes)))

#BP
#create topgo data object
myGOdata= new("topGOdata", description="My project", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO=geneID2GO)
#test for significance
#run weighted algorithm as classic doesnt take into consideration GO hierarchy so could overrepresent enrichment
resultFisher = runTest(myGOdata, algorithm="weight01", statistic="fisher")
#generate a table of results using Gentable function
allGO = usedGO(object = myGOdata)
allRes = GenTable(myGOdata, weightFisher = resultFisher, orderBy = "resultsFisher", ranksOf = "weightFisher", topNodes = length(allGO))
#change format to non-scientific digits
options(scipen = 999)
#correct for multiple testing e.g. a p-adjusted value
allRes$adjusted.p= p.adjust(allRes$weightFisher, method = "bonferroni", n = length(allRes$weightFisher))
allRes$q.value= p.adjust(allRes$weightFisher, method = "fdr", n = length(allRes$weightFisher))
BP=allRes
BP$ontology="BP"

#MF
myGOdata= new("topGOdata", description="My project", ontology="MF", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO=geneID2GO)
resultFisher = runTest(myGOdata, algorithm="weight01", statistic="fisher")
allGO = usedGO(object = myGOdata)
allRes = GenTable(myGOdata, weightFisher = resultFisher, orderBy = "resultFisher", ranksOf = "weightFisher", topNodes = length(allGO))
options(scipen = 999)
allRes$adjusted.p= p.adjust(allRes$weightFisher, method = "bonferroni", n = length(allRes$weightFisher))
allRes$q.value= p.adjust(allRes$weightFisher, method = "fdr", n = length(allRes$weightFisher))
MF=allRes
MF$ontology="MF"

#CC
myGOdata= new("topGOdata", description="My project", ontology="CC", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO=geneID2GO)
resultFisher = runTest(myGOdata, algorithm="weight01", statistic="fisher")
allGO = usedGO(object = myGOdata)
allRes = GenTable(myGOdata, weightFisher = resultFisher, orderBy = "resultFisher", ranksOf = "weightFisher", topNodes = length(allGO))
options(scipen = 999)
allRes$adjusted.p= p.adjust(allRes$weightFisher, method = "bonferroni", n = length(allRes$weightFisher))
allRes$q.value= p.adjust(allRes$weightFisher, method = "fdr", n = length(allRes$weightFisher))
CC=allRes
CC$ontology="CC"

out=rbind(BP, CC, MF)
out.s=subset(out, out$weightFisher <0.001)
write.table(out.s, "Results/GO_results_DE_FC0.0_ahya_34_21vs34_22.fisher_DONE.txt", row.names = FALSE, sep = "\t", quote = FALSE) #change accordingly

################### heatgrid GO terms ################################

#subset of interesting GO_enriched terms between heat+hypoxia vs heat only
GO_list=read.delim("DEG_GO_enriched_heatgrid.txt", sep="\t")

#manually curated categories of GO enriched terms
a<-c("Energy", "Immunity", "ROS", "Stress", "Iron")

GO_list$type<-factor(GO_list$type, levels=a)
b<-c("30_21vs30_22", "34_21vs34_22", "30_21vs34_21", "30_22vs34_22")
GO_list$comp<-factor(GO_list$comp, levels=b)

ggplot(GO_list, aes(x=comp, y=Term, fill=weightFisher)) + geom_tile(color="gray90", size=0.05)+
  guides(fill=guide_legend(title="P-value"))+
  scale_fill_gradientn(colours=c("turquoise4", "turquoise1", "paleturquoise1"))+
  theme_bw()+ 
  scale_y_discrete(position = 'right')+
  theme(axis.text.x = element_blank())+
  theme(legend.title = element_text(size=10), axis.title=element_blank())+
  facet_grid(type~comp, scales="free",space="free", switch="y", labeller=label_wrap_gen(width=20, multi_line = TRUE))+
  theme(strip.text.y = element_text(angle = 0), strip.background.x = element_rect(color="gray20"))+
  theme(strip.text.y.left = element_text(angle = 0))+
  theme(strip.text.x = element_text(angle = 45))+
  theme(axis.text.y=element_text(size=10, hjust=0))+
  theme(panel.spacing.x=unit(0, "lines"))+
  theme(legend.text.align = 0, legend.text = element_text(size=8))

#ggsave("Heatgrid_GO_plot.pdf", dpi=300,width = 6, height = 8)
#final label edits in affinitydesigner

################### Venn diagrams of DEGs ################################

#heat stress comparisons 
DE1=read.table("Results/DE_output_ahya_30_21v34_21_FC0.0_p0.05.txt", sep="\t", header=T, row.names = NULL)
g1=DE1$row.names
DE2=read.table("Results/DE_output_ahya_30_22v34_22_FC0.0_p0.05.txt", sep="\t", header=T, row.names = NULL)
g2=DE2$row.names

#baseline comparisons
DE3=read.table("Results/DE_output_ahya_30_21v30_22_FC0.0_p0.05.txt", sep="\t", header=T, row.names = NULL)
g3=DE3$row.names
DE4=read.table("Results/DE_output_ahya_34_21v34_22_FC0.0_p0.05.txt", sep="\t", header=T, row.names = NULL)
g4=DE4$row.names

x2<-list(g3, g4) #change accordingly

ggVennDiagram(x2, label_alpha = 0, lty=1, edge_size= 2, set_size= 2, category.names = c("30_21v34_21", "30_22v34_22")) + 
  scale_colour_manual(values= c("purple", "purple")) +
  ggplot2::scale_fill_gradient(low="white",high = "#CCCCFF")+ 
  theme(legend.position = "none")

#ggsave("venn_30_21v34_21with30_22v34_22.pdf", width=3, height=3, dpi = 300)
#ggsave("venn_30_21v30_22with34_21v34_22.pdf", width=3, height=3, dpi = 300)

################## annotate the overlapping genes #######################

#baseline donor (2021) vs nursery (2022) overlapping

#subset common DEGs
common_B <- subset(DE3, DE3$row.names %in% DE4$row.names)
#annotate DEGs based on eggnog annotations
anno_common_B<- subset(annot, annot$X.query %in% common_B$row.names)

#merge annotation with ddseq2 results
colnames(common_B)[which(names(common_B) == "row.names")] <- "Gene_name"
colnames(anno_common_B)[which(names(anno_common_B) == "X.query")] <- "Gene_name"
common_B_edited<-semi_join(common_B, anno_common_B, by = c("Gene_name"))
anno_common_B_dups <- anno_common_B[!duplicated(anno_common_B$Gene_name), ]

combined<-merge(x=anno_common_B_dups,y=common_B_edited, 
                by=c("Gene_name", "Gene_name"))

#subset genes of interest
histone2 <- subset(combined, combined$Description %like% "histone")
histone$Type = rep("Histone",length(nrow(histone)))
telomere <- subset(combined, combined$Description %like% "telo")
telomere$Type = rep("Telomere",length(nrow(telomere)))
chromatin <- subset(combined, combined$Description %like% "chromatin")
chromatin$Type = rep("Chromatin",length(nrow(chromatin)))

combined_GOI=rbind(histone, telomere, chromatin)

write.table(combined_GOI, "Results/Venn_overlap_DEG_30_21v30_22with34_21v34_22.txt", row.names = FALSE, sep = "\t", quote = FALSE) #change accordingly

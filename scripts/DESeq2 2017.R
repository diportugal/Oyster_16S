#Analysis of 2017 Data With DESeq2
#Preparing 2017 data for DESeq2
#2021-07-10
#Author:Diana Portugal
#Contact: dportugal8@gmail.com 



#Loading libraries ####
library("phyloseq"); packageVersion("phyloseq")
library("DESeq2")
packageVersion("DESeq2")
library("GenomicRanges")
library("GenomeInfoDb") 
require(DESeq2) #Not loading properly
#install.packages("DESeq2")#I think this fixed it....
library(dplyr)
library(tidyverse)



#Data re-load ####
OTU <- otu_table(otumat_matrix, taxa_are_rows = FALSE)
TAX <- tax_table(taxmat_matrix)
SAMP <- sample_data(metadata17_df)


physeq_class = phyloseq(OTU, TAX, SAMP) 
physeq_class


#Start here ####

#Analysis of RFTM Scores ####
#Make RFTM_score.x not a factor - Make it a numerical value 
##DESeq2 Plot of RFTM Score filtered by Phylum and Genus ####
RFTM_dds17 <- phyloseq_to_deseq2(physeq_class, ~RFTM_score.x) 
RFTM_dds17 <- DESeq(RFTM_dds17, test="Wald", fitType="parametric")

RFTM_res17= results(RFTM_dds17, cooksCutoff = FALSE)
alpha = 0.05
RFTM_sig17 = RFTM_res17[which(RFTM_res17$padj < alpha), ]
RFTM_sig17 = cbind(as(RFTM_sig17, "data.frame"), as(tax_table(physeq_class)[rownames(RFTM_sig17), ], "matrix"))

head(RFTM_sig17)
dim(RFTM_sig17)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)}

#Phylum = Variable 1
x = tapply(RFTM_sig17$log2FoldChange, RFTM_sig17$Phylum, function(x) max(x))
x = sort(x, TRUE)
RFTM_sig17$Phylum = factor(as.character(RFTM_sig17$Phylum), levels= names(x))

#Class = Variable 2
x = tapply(RFTM_sig17$log2FoldChange, RFTM_sig17$Genus.x, function(x) max(x))
x = sort(x, TRUE)
RFTM_sig17$Genus = factor(as.character(RFTM_sig17$Genus.x), levels= names(x))


ggplot(RFTM_sig17, aes(x=Genus.x, y=log2FoldChange, color=Phylum))+
  geom_point(size=2)+ 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5,))+
  labs(title = "DESeq2 of RFTM/OTU",
       subtitle = "Comparing Genus and Phylum presense as it compares to RFTM score",
       caption = "Data source: Oyster 16s 2017 (Using new taxa data)")+
  scale_colour_manual(values=c("#FA7169", "#33BEFF", "#8FC172","#FAC069", "#BA96D9", "#FCCFF4", "#BCB0EE", "#A6D5FD", "#4A456A", "#F7F3CD", "#F2B6AE" ))


##DESeq2 Plot of RFTM Score filtered by Phylum and Phylum ####

RFTM_dds17 <- phyloseq_to_deseq2(physeq_class, ~RFTM_score.x) 
RFTM_dds17 <- DESeq(RFTM_dds17, test="Wald", fitType="parametric")

RFTM_res17= results(RFTM_dds17, cooksCutoff = FALSE)
alpha = 0.05
RFTM_sig17 = RFTM_res17[which(RFTM_res17$padj < alpha), ]
RFTM_sig17 = cbind(as(RFTM_sig17, "data.frame"), as(tax_table(physeq_class)[rownames(RFTM_sig17), ], "matrix"))

head(RFTM_sig17)
dim(RFTM_sig17)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)}

#Phylum = Variable 1
x = tapply(RFTM_sig17$log2FoldChange, RFTM_sig17$Phylum, function(x) max(x))
x = sort(x, TRUE)
RFTM_sig17$Phylum = factor(as.character(RFTM_sig17$Phylum), levels= names(x))

#Class = Variable 2
x = tapply(RFTM_sig17$log2FoldChange, RFTM_sig17$Genus, function(x) max(x))
x = sort(x, TRUE)
RFTM_sig17$Genus = factor(as.character(RFTM_sig17$Genus), levels= names(x))


ggplot(RFTM_sig17, aes(x=Genus, y=log2FoldChange, color=Phylum))+geom_point(size=2)+ 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5,))+
  labs(title = "DESeq2 of RFTM/OTU",
       subtitle = "Comparing Phylum and Phylum presense as it compares to RFTM score",
       caption = "Data source: Oyster 16s 2017")+
  scale_colour_manual(values=c("#FA7169", "#33BEFF", "#8FC172","#FAC069", "#BA96D9", "#FCCFF4", "#BCB0EE", "#A6D5FD", "#4A456A", "#F7F3CD", "#F2B6AE" ))

ggsave("DESeq2_RFTM_Phylum_Genus.jpeg",width = 7, height = 5)



##DESeq2 Plot of RFTM Score filtered by Phylum and Class ####
ggplot(RFTM_sig17, aes(x=Class, y=log2FoldChange, color=Phylum))+
  geom_point(size=2)+ 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5,))+
  labs(title = "DESeq2 of RFTM/OTU",
       subtitle = "Comparing Class and Phylum presense as it compares to RFTM score",
       caption = "Data source: Oyster 16s 2017")+
  scale_colour_manual(values=c("#FA7169", "#33BEFF", "#8FC172","#FAC069", "#BA96D9", "#FCCFF4", "#BCB0EE", "#A6D5FD", "#4A456A", "#F7F3CD", "#F2B6AE" ))

ggsave("DESeq2_RFTM_Phylum_Class.jpeg",width = 7, height = 5)


##DESeq2 Plot of RFTM Score filtered by Phylum and Order ####
ggplot(RFTM_sig17, aes(x=Order, y=log2FoldChange, color=Phylum))+geom_point(size=2)+ 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5,))+
  labs(title = "DESeq2 of RFTM/OTU",
       subtitle = "Comparing Order and Phylum presense as it compares to RFTM score",
       caption = "Data source: Oyster 16s 2017")+
  scale_colour_manual(values=c("#FA7169", "#33BEFF", "#8FC172","#FAC069", "#BA96D9", "#FCCFF4", "#BCB0EE", "#A6D5FD", "#4A456A", "#F7F3CD", "#F2B6AE" ))

ggsave("DESeq2_RFTM_Phylum_Order.jpeg",width = 7, height = 5)


##DESeq2 Plot of RFTM Score filtered by Phylum and Family ####
ggplot(RFTM_sig17, aes(x=Family, y=log2FoldChange, color=Phylum))+geom_point(size=2)+ 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5,))+
  labs(title = "DESeq2 of RFTM/OTU",
       subtitle = "Comparing Family and Phylum presense as it compares to RFTM score",
       caption = "Data source: Oyster 16s 2017")+
  scale_colour_manual(values=c("#FA7169", "#33BEFF", "#8FC172","#FAC069", "#BA96D9", "#FCCFF4", "#BCB0EE", "#A6D5FD", "#4A456A", "#F7F3CD", "#F2B6AE" ))

ggsave("DESeq2_RFTM_Phylum_Family.jpeg",width = 7, height = 5)




#Pause here ####
##Analysis of Treatment2 #### 

Treat2_dds17 <- phyloseq_to_deseq2(physeq_class, ~Treatment2) 
Treat2_dds17 <- DESeq(Treat2_dds17, test="Wald", fitType="parametric")

Treat2_res17= results(Treat2_dds17, cooksCutoff = FALSE)
alpha = 0.05
Treat2_sig17 = Treat2_res17[which(Treat2_res17$padj < alpha), ]
Treat2_sig17 = cbind(as(Treat2_sig17, "data.frame"), as(tax_table(physeq_class)[rownames(Treat2_sig17), ], "matrix"))

head(Treat2_sig17)
dim(Treat2_sig17)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)}

#Phylum = Variable 1
x = tapply(Treat2_sig17$log2FoldChange, Treat2_sig17$Phylum, function(x) max(x))
x = sort(x, TRUE)
Treat2_sig17$Phylum = factor(as.character(Treat2_sig17$Phylum), levels= names(x))

#Genus = Variable 2
x = tapply(Treat2_sig17$log2FoldChange, Treat2_sig17$Genus, function(x) max(x))
x = sort(x, TRUE)
Treat2_sig17$Genus = factor(as.character(Treat2_sig17$Genus), levels= names(x))


Treat2_DESeq2 <- ggplot(Treat2_sig17, aes(x=Genus, y=log2FoldChange, color=Phylum))+geom_point(size=3)+ 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5,))+
  labs(title = "DESeq2 of Treatment2/OTU",
       subtitle = "Comparing phyla and genus presense as it compares to Treatment2",
       caption = "Data source: Oyster 16s 2017")+
  scale_colour_manual(values=c("#FA7169", "#33BEFF", "#8FC172","#FAC069", "#BA96D9", "#FCCFF4", "#BCB0EE", "#A6D5FD", "#4A456A", "#F7F3CD", "#F2B6AE" ))

print(Treat2_DESeq2)

ggsave("Treat2_DESeq2.jpeg", width = 7, height = 5)  






##Additional plot types ####

plotMA( res, ylim = c(-1, 1) ) #what is this?
#The x axis is the average expression over all samples, the y axis the log2 fold change between treatment and control.
#Genes with an adjusted p value below a threshold (here 0.1, the default) are shown in blue

plotDispEsts( rftm_dds, ylim = c(1e-6, 1e1) )

plotPCA( sigtab, intgroup = c( "Treatment2", "RFTM_score.x"), col=cols )
























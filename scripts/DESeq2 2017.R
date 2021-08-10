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
require(DESeq2) 
#install.packages("DESeq2")
library(dplyr)
library(tidyverse)
#install.packages("genefilter")
require(genefilter)



#Data re-load ####
OTU <- otu_table(otumat_matrix, taxa_are_rows = FALSE)
TAX <- tax_table(taxmat_matrix)
SAMP <- sample_data(metadata17_df)


physeq_class = phyloseq(OTU, TAX, SAMP) 
physeq_class


#Start here ####

#Analysis of RFTM Scores ####
#Make RFTM_score.x not a factor - Make it a numerical value 
RFTM_dds17 <- phyloseq_to_deseq2(physeq_class, ~RFTM_score.x) 
RFTM_dds17 <- DESeq(RFTM_dds17, test="Wald", fitType="parametric")

RFTM_res17= results(RFTM_dds17, cooksCutoff = FALSE)
alpha = 0.05
RFTM_sig17 = RFTM_res17[which(RFTM_res17$padj < alpha), ]
RFTM_sig17 = cbind(as(RFTM_sig17, "data.frame"), as(tax_table(physeq_class)[rownames(RFTM_sig17), ], "matrix"))

head(RFTM_sig17)
dim(RFTM_sig17)


#Phylum = Variable 1
x = tapply(RFTM_sig17$log2FoldChange, RFTM_sig17$Phylum, function(x) max(x))
x = sort(x, TRUE)
RFTM_sig17$Phylum = factor(as.character(RFTM_sig17$Phylum), levels= names(x))

#Class = Variable 2
x = tapply(RFTM_sig17$log2FoldChange, RFTM_sig17$Genus.x, function(x) max(x))
x = sort(x, TRUE)
RFTM_sig17$Genus = factor(as.character(RFTM_sig17$Genus.x), levels= names(x))


##Kingdom ####
ggplot(RFTM_sig17, aes(x=Kingdom, y=log2FoldChange, color=Phylum))+
  geom_point(size=2)+ 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5,))+
  labs(title = "DESeq2 of RFTM",
       subtitle = "Comparing Kingdom and Phylum presense as it compares to RFTM score",
       caption = "Data source: Oyster 16s 2017")
ggsave(filename = "DES_RFTM_Kingdom.jpeg", plot=last_plot(), path ="Plot Diagrams/DESeq2/", width = 7, height = 5)  


## Phylum ####
ggplot(RFTM_sig17, aes(x=Phylum, y=log2FoldChange, color=Phylum))+
  geom_point(size=2)+ 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5,))+
  labs(title = "DESeq2 of RFTM",
       subtitle = "Comparing Phylum and Phylum presense as it compares to RFTM score",
       caption = "Data source: Oyster 16s 2017")
ggsave(filename = "DES_RFTM_Phylum.jpeg", plot=last_plot(), path ="Plot Diagrams/DESeq2/", width = 7, height = 5)  


## Class ####
ggplot(RFTM_sig17, aes(x=Class, y=log2FoldChange, color=Phylum))+
  geom_point(size=2)+ 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5,))+
  labs(title = "DESeq2 of RFTM",
       subtitle = "Comparing Class and Phylum presense as it compares to RFTM score",
       caption = "Data source: Oyster 16s 2017")
ggsave(filename = "DES_RFTM_Class.jpeg", plot=last_plot(), path ="Plot Diagrams/DESeq2/", width = 7, height = 5)  


##Order ####
ggplot(RFTM_sig17, aes(x=Order, y=log2FoldChange, color=Phylum))+
  geom_point(size=2)+ 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5,))+
  labs(title = "DESeq2 of RFTM",
       subtitle = "Comparing Order and Phylum presense as it compares to RFTM score",
       caption = "Data source: Oyster 16s 2017")
ggsave(filename = "DES_RFTM_Order.jpeg", plot=last_plot(), path ="Plot Diagrams/DESeq2/", width = 7, height = 5)  



##Family ####
ggplot(RFTM_sig17, aes(x=Family, y=log2FoldChange, color=Phylum))+
  geom_point(size=2)+ 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5,))+
  labs(title = "DESeq2 of RFTM",
       subtitle = "Comparing Family and Phylum presense as it compares to RFTM score",
       caption = "Data source: Oyster 16s 2017")
ggsave(filename = "DES_RFTM_Family.jpeg", plot=last_plot(), path ="Plot Diagrams/DESeq2/", width = 7, height = 5)  


##Genus ####
ggplot(RFTM_sig17, aes(x=Genus.x, y=log2FoldChange, color=Phylum))+
  geom_point(size=2)+ 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5,))+
  labs(title = "DESeq2 of RFTM",
       subtitle = "Comparing Genus and Phylum presense as it compares to RFTM score",
       caption = "Data source: Oyster 16s 2017")
ggsave(filename = "DES_RFTM_Genus.jpeg", plot=last_plot(), path ="Plot Diagrams/DESeq2/", width = 7, height = 5)  



#*####
#DESEQ2 of PeaCrabs ####
#Analysis of Peacrab Scores ####

PEA_dds17 <- phyloseq_to_deseq2(physeq_class, ~peacrabs.x) 
PEA_dds17 <- DESeq(PEA_dds17, test="Wald", fitType="parametric")

PEA_res17= results(PEA_dds17, cooksCutoff = FALSE)
alpha = 0.05
PEA_sig17 = PEA_res17[which(PEA_res17$padj < alpha), ]
PEA_sig17 = cbind(as(PEA_sig17, "data.frame"), as(tax_table(physeq_class)[rownames(PEA_sig17), ], "matrix"))

head(PEA_sig17)
dim(PEA_sig17)

#Phylum = Variable 1
x = tapply(PEA_sig17$log2FoldChange, PEA_sig17$Phylum, function(x) max(x))
x = sort(x, TRUE)
PEA_sig17$Phylum = factor(as.character(PEA_sig17$Phylum), levels= names(x))

#Genus = Variable 2
x = tapply(PEA_sig17$log2FoldChange, PEA_sig17$Family, function(x) max(x))
x = sort(x, TRUE)
PEA_sig17$Family = factor(as.character(PEA_sig17$Family), levels= names(x))


##Kingdom ####
ggplot(PEA_sig17, aes(x=Kingdom, y=log2FoldChange, color=Phylum))+
  geom_point(size=2)+ 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5,))+
  labs(title = "DESeq2 of Peacrab",
       subtitle = "Comparing Kingdom and Phylum presense as it compares to Peacrab prescense",
       caption = "Data source: Oyster 16s 2017")
ggsave("DES_PEA_Kingdom.jpeg", path="Plot Diagrams/DESeq2/", width = 7, height = 5)


##Phylum ####
ggplot(PEA_sig17, aes(x=Phylum, y=log2FoldChange, color=Phylum))+
  geom_point(size=2)+ 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5,))+
  labs(title = "DESeq2 of Peacrab",
       subtitle = "Comparing Phylum and Phylum presense as it compares to Peacrab prescense",
       caption = "Data source: Oyster 16s 2017")
ggsave("DES_PEA_Phylum.jpeg", path="Plot Diagrams/DESeq2/", width = 7, height = 5)


##Class ####
ggplot(PEA_sig17, aes(x=Class, y=log2FoldChange, color=Phylum))+
  geom_point(size=2)+ 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5,))+
  labs(title = "DESeq2 of Peacrab",
       subtitle = "Comparing Class and Phylum presense as it compares to Peacrab prescense",
       caption = "Data source: Oyster 16s 2017")
ggsave("DES_PEA_Class.jpeg", path="Plot Diagrams/DESeq2/", width = 7, height = 5)


##Order ####
ggplot(PEA_sig17, aes(x=Order, y=log2FoldChange, color=Phylum))+
  geom_point(size=2)+ 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5,))+
  labs(title = "DESeq2 of Peacrab",
       subtitle = "Comparing Order and Phylum presense as it compares to Peacrab prescense",
       caption = "Data source: Oyster 16s 2017")
ggsave("DES_PEA_Order.jpeg", path="Plot Diagrams/DESeq2/", width = 7, height = 5)


##Family ####
ggplot(PEA_sig17, aes(x=Family, y=log2FoldChange, color=Phylum))+
  geom_point(size=2)+ 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5,))+
  labs(title = "DESeq2 of Peacrab",
       subtitle = "Comparing Family and Phylum presense as it compares to Peacrab prescense",
       caption = "Data source: Oyster 16s 2017")
ggsave("DES_PEA_Family.jpeg", path="Plot Diagrams/DESeq2/", width = 7, height = 5)


##Genus ####
ggplot(PEA_sig17, aes(x=Genus.x, y=log2FoldChange, color=Phylum))+
  geom_point(size=2)+ 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5,))+
  labs(title = "DESeq2 of Peacrab",
       subtitle = "Comparing Genus and Phylum presense as it compares to Peacrab prescense",
       caption = "Data source: Oyster 16s 2017")
ggsave("DES_PEA_Genus.jpeg", path="Plot Diagrams/DESeq2/", width = 7, height = 5)

#*####



#Additional plot types ####

plotMA( res, ylim = c(-1, 1) ) 
#The x axis is the average expression over all samples, the y axis the log2 fold change between treatment and control.
#Genes with an adjusted p value below a threshold (here 0.1, the default) are shown in blue

plotDispEsts( rftm_dds, ylim = c(1e-6, 1e1) )

plotPCA( sigtab, intgroup = c( "Treatment2", "RFTM_score.x"), col=cols )





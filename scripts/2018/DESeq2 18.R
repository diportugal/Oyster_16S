#Analysis of 2018 Data With DESeq2
#Preparing 2018 data for DESeq2
#2021-07-14
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
OTU18 <- otu_table(otumat_matrix18, taxa_are_rows = FALSE)
TAX <- tax_table(taxmat_matrix)
SAMP18 <- sample_data(metadata18_df)


physeq_class18 = phyloseq(OTU18, TAX, SAMP18) 
physeq_class18



#Setting up the data ####

#Start here ####

#Analysis of RFTM Scores ####
#Make RFTM_score.x not a factor - Make it a numerical value 
#dds = estimateSizeFactors(dds, geoMeans=geoMeans, locfunc=shorth

gm_mean <-function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0])))
geoMeans <- apply(OTU18, 2, gm_mean)
RFTM_dds18 <- phyloseq_to_deseq2(physeq_class18, ~RFTM_score.x) 
RFTM_dds18 = estimateSizeFactors(RFTM_dds18, geoMeans=geoMeans, locfunc=shorth)

view(OTU18)
dim(OTU18)
dim(TAX)
dim(RFTM_dds18)
nrow(RFTM_dds18)
dim(geoMeans)
view(geoMeans)



RFTM_dds18 <- DESeq(RFTM_dds18, test="Wald", fitType="parametric")

RFTM_res18= results(RFTM_dds18, cooksCutoff = FALSE)
alpha = 0.05
RFTM_sig18 = RFTM_res18[which(RFTM_res18$padj < alpha), ]
RFTM_sig18 = cbind(as(RFTM_sig18, "data.frame"), as(tax_table(physeq_class18)[rownames(RFTM_sig18), ], "matrix"))

head(RFTM_sig18)
dim(RFTM_sig18)


#Phylum = Variable 1
x = tapply(RFTM_sig18$log2FoldChange, RFTM_sig18$Phylum, function(x) max(x))
x = sort(x, TRUE)
RFTM_sig18$Phylum = factor(as.character(RFTM_sig18$Phylum), levels= names(x))

#Class = Variable 2
x = tapply(RFTM_sig18$log2FoldChange, RFTM_sig18$Class, function(x) max(x))
x = sort(x, TRUE)
RFTM_sig18$Class = factor(as.character(RFTM_sig18$Class), levels= names(x))


ggplot(RFTM_sig18, aes(x=Class, y=log2FoldChange, color=Phylum))+
  geom_point(size=3)+ 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5,))+
  labs(title = "DESeq2 of RFTM/OTU",
       subtitle = "Comparing Phylum and Class presense as it compares to RFTM score",
       caption = "Data source: Oyster 16s 2018")

ggsave("DESeq2_RFTM_Phylum_Class.jpeg",width = 7, height = 5)



##DESeq2 Plot of RFTM Score filtered by Phylum and Order ####
ggplot(RFTM_sig18, aes(x=Order, y=log2FoldChange, color=Phylum))+
  geom_point(size=3)+ 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5,))+
  labs(title = "DESeq2 of RFTM/OTU",
       subtitle = "Comparing Phylum and Order presense as it compares to RFTM score",
       caption = "Data source: Oyster 16s 2018")

ggsave("DESeq2_RFTM_Phylum_Order.jpeg",width = 7, height = 5)



##DESeq2 Plot of RFTM Score filtered by Phylum and Family ####
ggplot(RFTM_sig18, aes(x=Family, y=log2FoldChange, color=Phylum))+
  geom_point(size=3)+ 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5,))+
  labs(title = "DESeq2 of RFTM/OTU",
       subtitle = "Comparing Phylum and Family presense as it compares to RFTM score",
       caption = "Data source: Oyster 16s 2018")

ggsave("DESeq2_RFTM_Phylum_Family.jpeg",width = 7, height = 5)



##DESeq2 Plot of RFTM Score filtered by Phylum and Genus.x ####
ggplot(RFTM_sig18, aes(x=Genus.x, y=log2FoldChange, color=Phylum))+
  geom_point(size=3)+ 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5,))+
  labs(title = "DESeq2 of RFTM/OTU",
       subtitle = "Comparing Phylum and Genus.x presense as it compares to RFTM score",
       caption = "Data source: Oyster 16s 2018")

ggsave("DESeq2_RFTM_Phylum_Genus.x.jpeg",width = 7, height = 5)









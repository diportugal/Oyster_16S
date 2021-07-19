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
require(DESeq2) #Not loading properly
install.packages("DESeq2")#I think this fixed it....
library(dplyr)
library(tidyverse)


#Data re-load ####
OTU18 <- otu_table(otumat_matrix18, taxa_are_rows = FALSE)
TAX <- tax_table(taxmat_matrix)
SAMP18 <- sample_data(metadata18_df)


physeq_class18 = phyloseq(OTU18, TAX, SAMP18) 
physeq_class18



#Setting up the data ####

physeq_class18


#Start here ####

#Analysis of RFTM Scores ####
#Make RFTM_score.x not a factor - Make it a numerical value 
#dds = estimateSizeFactors(dds, geoMeans=geoMeans, locfunc=shorth


RFTM_dds18 <- phyloseq_to_deseq2(physeq_class18, ~RFTM_score.x) 
RFTM_dds18 <- estimateSizeFactors(RFTM_dds18, locfunc<-"shorth")
RFTM_dds18 <- DESeq(RFTM_dds18, test="Wald", fitType="parametric")

RFTM_res18= results(RFTM_dds18, cooksCutoff = FALSE)
alpha = 0.05
RFTM_sig18 = RFTM_res18[which(RFTM_res18$padj < alpha), ]
RFTM_sig18 = cbind(as(RFTM_sig18, "data.frame"), as(tax_table(physeq_class18)[rownames(RFTM_sig18), ], "matrix"))

head(RFTM_sig18)
dim(RFTM_sig18)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)}

#Phylum = Variable 1
x = tapply(RFTM_sig18$log2FoldChange, RFTM_sig18$Phylum, function(x) max(x))
x = sort(x, TRUE)
RFTM_sig18$Phylum = factor(as.character(RFTM_sig18$Phylum), levels= names(x))

#Genus = Variable 2
x = tapply(RFTM_sig18$log2FoldChange, RFTM_sig18$Genus, function(x) max(x))
x = sort(x, TRUE)
RFTM_sig18$Genus = factor(as.character(RFTM_sig18$Genus), levels= names(x))


RFTM_DESeq2_18 <- ggplot(RFTM_sig18, aes(x=Genus, y=log2FoldChange, color=Phylum))+geom_point(size=3)+ 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5,))+
  labs(title = "DESeq2 of RFTM/OTU",
       subtitle = "Comparing phyla and genus presense as it compares to RFTM score",
       caption = "Data source: Oyster 16s 2018")+
  scale_colour_manual(values=c("#FA7169", "#33BEFF", "#8FC172","#FAC069", "#BA96D9", "#FCCFF4", "#BCB0EE", "#A6D5FD", "#4A456A", "#F7F3CD", "#F2B6AE" ))

#ggsave("RFTM_DESeq2.jpeg", width = 7, height = 5)  

























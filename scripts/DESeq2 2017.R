# Analysis of 2017 Data With DESeq2
#Preparing 2017 data for DESeq2
#2021-07-10
#Author:Diana Portugal
#Contact: dportugal8@gmail.com 



#Loading libraries ####
library("phyloseq"); packageVersion("phyloseq")
library("DESeq2")
packageVersion("DESeq2")
library("GenomicRanges")
library("GenomeInfoDb") #ERROR LOADING ####
require(DESeq2) #Not loading properly
install.packages("DESeq2")#I think this fixed it....
library(dplyr)
library(tidyverse)



#Data re-load ####
OTU <- otu_table(otumat_matrix, taxa_are_rows = FALSE)
TAX <- tax_table(taxmat_matrix)
SAMP <- sample_data(metadata17_df)


physeq_class = phyloseq(OTU, TAX, SAMP) 
physeq_class


#Start here ####
#Simple plot ####
physeq_class
head(sample_data(physeq_class)$Treatment2)
physeq_class= subset_samples(physeq_class, )


deseq.17.Rscore = phyloseq_to_deseq2(physeq_class, ~ Site.x)
deseq.17.Rscore = DESeq(deseq.17.Rscore, test="Wald", fitType="parametric")


results.17.Rscore = results(deseq.17.Rscore, cooksCutoff = FALSE)
alpha = 0.05
sigtab.17.Rscore = results.17.Rscore[which(results.17.Rscore$padj < alpha), ]
sigtab.17.Rscore = cbind(as(sigtab.17.Rscore, "data.frame"), as(tax_table(physeq_class)[rownames(sigtab.17.Rscore), ], "matrix"))


theme_set(theme_bw())

x = tapply(sigtab.17.Rscore$log2FoldChange, sigtab.17.Rscore$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab.17.Rscore$Genus = factor(as.character(sigtab.17.Rscore$Genus), levels=names(x))


sigtab.17.Rscore$Class = factor(as.character(sigtab.17.Rscore$Genus), levels=names(x))

ggplot(sigtab.17.Rscore, aes(x=Order, y=log2FoldChange, color= Genus)) + geom_point(size=6)



#Test 2 ####
m <- phyloseq_to_deseq2(physeq_class, ~ RFTM_score.x)
m <- DESeq(m, test="Wald", fitType="parametric")


res.m.17 = results(m, cooksCutoff = FALSE)
alpha = 0.05
sigtab.m.17 = res.m.17[which(res.m.17$padj < alpha), ]
sigtab.m.17 = cbind(as(sigtab.m.17, "data.frame"), as(tax_table(physeq_class)[rownames(sigtab.m.17), ], "matrix"))

dim(sigtab.m.17)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
#Order Order
x = tapply(sigtab.m.17$log2FoldChange, sigtab.m.17$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab.m.17$Order = factor(as.character(sigtab.m.17$Order), levels= names(x))

#Genus Order
x = tapply(sigtab.m.17$log2FoldChange, sigtab.m.17$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab.m.17$Genus = factor(as.character(sigtab.m.17$Genus), levels= names(x))


ggplot(sigtab.m.17, aes(x=Genus, y=log2FoldChange, color=Order)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))



#Test 3 #### 

m <- phyloseq_to_deseq2(physeq_class, ~ RFTM_score.x)
m <- DESeq(m, test="Wald", fitType="parametric")


res.m.17 = results(m, cooksCutoff = FALSE)
alpha = 0.05
sigtab.m.17 = res.m.17[which(res.m.17$padj < alpha), ]
sigtab.m.17 = cbind(as(sigtab.m.17, "data.frame"), as(tax_table(physeq_class)[rownames(sigtab.m.17), ], "matrix"))

dim(sigtab.m.17)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
#Order Order
x = tapply(sigtab.m.17$log2FoldChange, sigtab.m.17$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab.m.17$Order = factor(as.character(sigtab.m.17$Order), levels= names(x))

#Genus Order
x = tapply(sigtab.m.17$log2FoldChange, sigtab.m.17$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab.m.17$Genus = factor(as.character(sigtab.m.17$Genus), levels= names(x))


ggplot(sigtab.m.17, aes(x=Genus, y=log2FoldChange, color=Order)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))








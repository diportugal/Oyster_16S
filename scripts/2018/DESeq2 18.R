# Analysis of 2018 Data With DESeq2
#Preparing 2018 data for DESeq2
#2021-07-10
#Author:Diana Portugal
#Contact: dportugal8@gmail.com 



#Loading Libraries ####

library("phyloseq"); packageVersion("phyloseq")
BiocManager::install("DESeq2")

library("DESeq2")
packageVersion("DESeq2")


#Setting up the data ####

desq18 = phyloseq_to_deseq2(physeq_class18, ~Species.x)

desq18 = DESeq(desq18, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(kostic)[rownames(sigtab), ], "matrix"))
head(sigtab)



library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))


OTU18 





















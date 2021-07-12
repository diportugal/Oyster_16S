# Analysis of 2017 Data With DESeq2
#Preparing 2017 data for DESeq2
#2021-07-10
#Author:Diana Portugal
#Contact: dportugal8@gmail.com 



#Loading libraries ####
library("phyloseq"); packageVersion("phyloseq")
library("DESeq2")
packageVersion("DESeq2")


#Data re-load ####
OTU <- otu_table(otumat_matrix, taxa_are_rows = FALSE)
TAX <- tax_table(taxmat_matrix)
SAMP <- sample_data(metadata17_df)


physeq_class = phyloseq(OTU, TAX, SAMP) 
physeq_class


#Start here ####

physeq_class
head(sample_data(physeq_class)$Species.x)


deseq.17 = phyloseq_to_deseq2(physeq_class, ~ RFTM_score.x)
deseq.17 = DESeq(deseq.17, test="Wald", fitType="parametric")


RFTM_score.x = results(deseq.17, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(results17$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(kostic)[rownames(sigtab), ], "matrix"))
head(sigtab)


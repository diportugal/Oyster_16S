#Document Information ####
#Oyster 16S 2017 Sample Data Analysis 
#Author: Diana Portugal 
#Contact: dportugal8@gmail.com 


library(data.table)
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library(plyr)
library(dplyr)
theme_set(theme_bw())


#Loading the data 
#ASV as OTU
otumat <- fread("Oyster_data_raw/asvtable_de17.csv")
rownames(otumat) <- otumat$V1
otumat$V1 = NULL
otumat <- as.matrix(otumat)
class(otumat)

#Run23 as TAX
taxmat <- fread("Oyster_data_raw/Run23_taxa.csv")
rownames(taxmat) <- taxmat$V1
taxmat$V1 = NULL
colnames(taxmat) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
taxmat <- as.matrix(taxmat)
class(taxmat)


#Cleanmetadata17 as SAMP
sampmat <- fread("Oyster_data_raw/cleanmetadata17")
rownames(sampmat) <- sampmat$V1
sampmat$V1 = NULL
class(sampmat)


OTU <- otu_table(otumat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
OTU
TAX

physeq = phyloseq(OTU, TAX)
physeq
# This is the last thing that works, when taxa_are_rows = TRUE, it doesnt work when = FALSE
# Does not let me add the SAMP value


plot_bar(physeq, fill = "Family") #Does not come up with a plot 




OTU=transform_sample_counts(OTU, function(x) 1E6 * x/sum(x))




physeq <- phyloseq(OTU, TAX)
physeq 




#PLOT ORDINANCE ####







#OTU <- transform_sample_counts(OTU, function(x) 1E6 * x/sum(x))


#PLOT BAR ####

data("physeq")
gp.ch <- subset_taxa(GlobalPatterns, Phylum == "Spirochaetota")


















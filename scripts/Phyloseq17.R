# Analysis of 2017 Data With Phyloseq
#Preparing 2017 data for Phyloseq
#2021-06-21
#Author:Diana Portugal
#Contact: dportugal8@gmail.com 


#PREP WORK ####
## LOADING REQUIRED R PACKAGES ####
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("plyr"); packageVersion("plyr")
library("data.table")
theme_set(theme_bw())



## IMPORTING DATA ####
metadata17 <- read.csv("Oyster_data_raw/meta17cleaned")

asvtable17 <- fread("Oyster_data_raw/asvtable_de17.csv")

run23 <- read.csv("Oyster_data_raw/Run123_taxa_complete.csv")


## CHANGING ROW NAMES FOR EACH DATA SET ####
rownames(metadata17) = metadata17$X.1
metadata17$X.1=NULL
rownames(metadata17)
#ROW NAMES ARE THE UNIQUE IDs 

rownames(asvtable17) = asvtable17$V1
rownames(asvtable17)
asvtable17$V1= NULL
#ROW NAMES ARE THE UNIQUE IDs 

rownames(run23) = run23$Row.names
run23$Row.names = NULL        
rownames(run23)
#ROW NAMES ARE THE SEQUENCE 


## CONVERTING TO MATRICIES ####
metadata17_df <- as.data.frame(metadata17, rownames("X.1"))
rownames(metadata17_df)
#STILL UNIQUE ID

otumat_matrix <- as.matrix(asvtable17, rownames=rownames(asvtable17))
rownames(otumat_matrix)
#STILL UNIQUE ID 


rownames(run23)

taxmat_matrix <- as.matrix(run23) 
colnames(taxmat_matrix) <- c("X", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus.x","Genus.y","Species")
rownames(taxmat_matrix)
#STILL SEQUENCE 



## SETTING OTU, TAX, SAMP ####
OTU <- otu_table(otumat_matrix, taxa_are_rows = FALSE)

TAX <- tax_table(taxmat_matrix)

SAMP <- sample_data(metadata17_df)



## INSPECTING SAMPLE NAMES####
sample_names(SAMP)
sample_names(OTU)
sample_names(TAX)


## EVENING OUT THE DATA ####
OTU=transform_sample_counts(OTU, function(x) 1E6 * x/sum(x))


## READING THROUGH PHYLOSEQ ####
physeq_class = phyloseq(OTU, TAX, SAMP) 
physeq_class


#Final Object ####
physeq_class

#Comment









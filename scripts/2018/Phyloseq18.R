# Analysis of 2018 Data With Phyloseq
#Preparing 2018 data for Phyloseq
#2021-06-27
#Author:Diana Portugal
#Contact: dportugal8@gmail.com 


## LOADING REQUIRED R PACKAGES ####
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("plyr"); packageVersion("plyr")
library("data.table")
theme_set(theme_bw())




## IMPORTING DATA ####
metadata18 <- read.csv("Oyster_data_raw/cleanmetadata18")
metadata18

asvtable18 <- fread("Oyster_data_raw/asvtable_de18.csv")

run23 <- read.csv("Oyster_data_raw/Run23_taxa.csv")

## CHANGING ROW NAMES FOR EACH DATA SET ####
rownames(metadata18) = metadata18$UniqueID
metadata18$UniqueID=NULL
rownames(metadata18)
#ROW NAMES ARE THE UNIQUE IDs 

rownames(asvtable18) = asvtable18$V1
rownames(asvtable18)
asvtable18$V1= NULL
#ROW NAMES ARE THE UNIQUE IDs 

rownames(run23) = run23$X
run23$X = NULL        
rownames(run23)
#ROW NAMES ARE THE SEQUENCE 





## CONVERTING TO MATRICIES ####
metadata18_df <- as.data.frame(metadata18, rownames("UniqueID"))
rownames(metadata18_df)
#STILL UNIQUE ID

otumat_matrix18 <- as.matrix(asvtable18, rownames=rownames(asvtable18))
rownames(otumat_matrix18)
#STILL UNIQUE ID 

taxmat_matrix <- as.matrix(run23) 
colnames(taxmat_matrix) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
rownames(taxmat_matrix)
#STILL SEQUENCE 



## SETTING OTU, TAX, SAMP ####
OTU18 <- otu_table(otumat_matrix18, taxa_are_rows = FALSE)

TAX <- tax_table(taxmat_matrix)

SAMP18 <- sample_data(metadata18_df)



## INSPECTING SAMPLE NAMES####
sample_names(SAMP18)
sample_names(OTU18)
sample_names(TAX)




## EVENING OUT THE DATA ####
OTU18=transform_sample_counts(OTU18, function(x) 1E6 * x/sum(x))




## READING THROUGH PHYLOSEQ ####
physeq_class18 = phyloseq(OTU18, TAX, SAMP18) 
physeq_class18



sample_data(physeq_class18)$RFTM_score.x=as.factor(sample_data(physeq_class18)$RFTM_score.x)

sample_data(physeq_class18)$Weight=as.factor(sample_data(physeq_class18)$Weight)



#End here with the phyloseq prep work ####
















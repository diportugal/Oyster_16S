library("phyloseq"); packageVersion("phyloseq")

library("ggplot2"); packageVersion("ggplot2")
library("plyr"); packageVersion("plyr")

theme_set(theme_bw())


#IMPORTING DATA ####

metadata17 <- read.csv("Oyster_data_raw/meta17cleaned")

asvtable17 <- fread("Oyster_data_raw/asvtable_de17.csv")

run23 <- read.csv("Oyster_data_raw/Run23_taxa.csv")


#CHANGING ROW NAMES ####

rownames(metadata17) = metadata17$X.1
metadata17$X.1=NULL
rownames(metadata17)
#ROW NAMES ARE THE UNIQUE IDs 

rownames(asvtable17) = asvtable17$V1
rownames(asvtable17)
asvtable17$V1= NULL
#ROW NAMES ARE THE UNIQUE IDs 

rownames(run23) = run23$X
run23$X = NULL        
rownames(run23)
#ROW NAMES ARE THE SEQUENCE 



#SETTING TAXMAT & OTUMAT
asvtable17
rownames(asvtable17)

taxmat <- run23
taxmat
rownames(taxmat)


#CONVERT TO MATRIX ####
metadata17_df <- as.data.frame(metadata17, rownames("X.1"))
rownames(metadata17_df)

otumat_matrix <- as.matrix(asvtable17, rownames=rownames(asvtable17))
rownames(otumat_matrix)
#STILL UNIQUE ID 

taxmat_matrix <- as.matrix(taxmat) 
colnames(taxmat_matrix) <- c("Kingdom", "Phylum", "Class", "Order", "Family", 
                             "Genus")
rownames(taxmat_matrix)
#STILL SEQUENCE 


# SETTING OTU, TAX, SAMP ####
OTU <- otu_table(otumat_matrix, taxa_are_rows = FALSE)

TAX <- tax_table(taxmat_matrix)

SAMP <- sample_data(metadata17_df)


#INSPECTING SAMPLE NAMES####
sample_names(SAMP)
sample_names(OTU)
sample_names(TAX)


#EVENING OUT THE DATA ####
OTU=transform_sample_counts(OTU, function(x) 1E6 * x/sum(x))


#READING THROUGH PHYLOSEQ ####
physeq_class = phyloseq(OTU, TAX, SAMP) 
physeq_class

#Comment 
#Comment
#com



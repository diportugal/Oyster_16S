library("phyloseq"); packageVersion("phyloseq")

data <- read.csv("Oyster_data_raw/cleanmetadata17")
library("ggplot2"); packageVersion("ggplot2")
library("plyr"); packageVersion("plyr")

theme_set(theme_bw())


#IMPORTING DATA ####

metadata17 <- fread("Oyster_data_raw/cleanmetadata17")

asvtable17 <- fread("Oyster_data_raw/asvtable_de17.csv")

run23 <- fread("Oyster_data_raw/Run23_taxa.csv")


#CHANGING ROW NAMES ####

rownames(metadata17) <- metadata17$UniqueID 
metadata17$UniqueID=NULL
metadata17$X=NULL
rownames(metadata17)


rownames(asvtable17) <- asvtable17$V1 
rownames(asvtable17)


rownames(run23) <- run23$V1 
rownames(run23)



#SETTING TAXMAT & OTUMAT

otumat <- asvtable17

taxmat <- run23


#CONVERT TO MATRIX ####

metadata17_df <- as.data.frame(metadata17)

otumat_matrix <- as.matrix(otumat, rownames = "V1")

taxmat_matrix <- as.matrix(taxmat, rownames = "V1")
colnames(taxmat_matrix) <- c("Kingdom", "Phylum", "Class", "Order", "Family", 
                          "Genus")

# SETTING OTU, TAX, SAMP ####


OTU <- otu_table(otumat_matrix, taxa_are_rows = FALSE)

TAX <- tax_table(taxmat_matrix)

SAMP <- sample_data(metadata17_df)


OTU=transform_sample_counts(OTU, function(x) 1E6 * x/sum(x))


physeq_class = phyloseq(OTU, TAX, SAMP) #ERROR ####
physeq_class


#Is there an error in my cleanmetadata17 file??



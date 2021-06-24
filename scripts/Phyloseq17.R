# Analysis of 2017 Data With Phyloseq
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

run23 <- read.csv("Oyster_data_raw/Run23_taxa.csv")



## CHANGING ROW NAMES FOR EACH DATA SET ####
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



## CONVERTING TO MATRICIES ####
metadata17_df <- as.data.frame(metadata17, rownames("X.1"))
rownames(metadata17_df)
#STILL UNIQUE ID

otumat_matrix <- as.matrix(asvtable17, rownames=rownames(asvtable17))
rownames(otumat_matrix)
#STILL UNIQUE ID 

taxmat_matrix <- as.matrix(run23) 
colnames(taxmat_matrix) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
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




#STARTING DATA ANALYSIS ####

## MAKING A PLOT GRAPH ####
plot_bar(physeq_class, x= "Treatment2", y="RFTM_score.x", fill = "Species.x")
#NOTE - The graphs take a while to load 



## MAKING AN NMDS PLOT (SCATTED DOTS): JUST OTUs ####
Phy.ord <- ordinate(physeq_class, "NMDS", "bray")
plot1 = plot_ordination(physeq_class, Phy.ord, type="taxa", color="Kingdom", title="TAXA")
print(plot1)

Phy.ord <- ordinate(physeq_class, "NMDS", "bray")
plot1 = plot_ordination(physeq_class, Phy.ord, type="taxa", color="Kingdom", title="TAXA")
print(plot1)

#Grouping by Phylum and doing 3 plots stacked on top of each other
plot1 + facet_wrap(~Kingdom, 3) 

#QUESTION: CAN WE CHANGE THE COLORS?!?!?



## MAKING A PLOT: JUST SAMPLES ####
plot2 = plot_ordination(physeq_class, Phy.ord, type="samples", color= sample_variables(physeq), shape=NULL) 
plot2 + geom_polygon(aes(fill=NULL)) + geom_point(size=2, color= "red") + ggtitle("PLOT 2: JUST SAMPLES")
#LOOK MORE INTO MANIPULATING THESE VARIABLES



## MAKING A BIPLOT GRAPHIC ####
plot3 = plot_ordination(physeq_class, Phy.ord, type="biplot", color=NULL, shape=NULL, title="PLOT 3: BIPLOT GRAPHIC")
# Some stuff to modify the automatic shape scale, not nessary to include but look into!!
GP1.shape.names = get_taxa_unique(GP1, "Phylum")
GP1.shape <- 15:(15 + length(GP1.shape.names) - 1)
names(GP1.shape) <- GP1.shape.names
GP1.shape["samples"] <- 16
p3 + scale_shape_manual(values=GP1.shape)

print(plot3)


## MAKING A SPLIT GRAPHIC ####

plot4 = plot_ordination(physeq_class, Phy.ord, type="split", color=NULL, shape=NULL, label=NULL, title="PLOT 4 : SPLIT GRAPHIC") 
print(plot4)




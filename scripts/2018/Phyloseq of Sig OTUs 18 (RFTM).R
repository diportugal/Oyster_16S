#Oyster 16S - Phyloseq Anaysis for 2018 RFTM Significant Data
#Author: Diana Portugal 
#Contact: dportugal8@gmail.com 
#Date 21-7-2021

library("dplyr")
library("ggpubr")
library("phyloseq"); packageVersion("phyloseq")
data("GlobalPatterns")
library("ggplot2"); packageVersion("ggplot2")
library("plyr"); packageVersion("plyr")
library("data.table")
library("vegan")
#install.packages("ggpubr")


#PREP WORK ####

##Isolating significant OTUs for RFTM Score####
phylo_rftm18 <- RFTM_sig18 %>% 
  select(Kingdom, Phylum, Class, Order, Family, Genus.x)


metadata18 <- read.csv("Oyster_data_raw/cleanmetadata18")
metadata18

asvtable18 <- fread("Oyster_data_raw/asvtable_de18.csv")

run23 <- read.csv("Oyster_data_raw/Run123_taxa_complete.csv")


## CHANGING ROW NAMES FOR EACH DATA SET ####
rownames(metadata18) = metadata18$UniqueID
metadata18$UniqueID=NULL
rownames(metadata18)
#ROW NAMES ARE THE UNIQUE IDs 

rownames(asvtable18) = asvtable18$V1
rownames(asvtable18)
asvtable18$V1= NULL
#ROW NAMES ARE THE UNIQUE IDs 

rownames(run23) = run23$Row.names
run23$Row.names = NULL        
rownames(run23)
#ROW NAMES ARE THE SEQUENCE 


## CONVERTING TO MATRICIES ####
metadata18_df <- as.data.frame(metadata18, rownames("UniqueID"))
rownames(metadata18_df)
#STILL UNIQUE ID

otumat_matrix18 <- as.matrix(asvtable18, rownames=rownames(asvtable18))
rownames(otumat_matrix18)
#STILL UNIQUE ID 

#Using significant otu data from DESeq2 analysis of RFTM score
sig_OTUs_RFTM18 <- as.matrix(phylo_rftm18)
colnames(phylo_rftm18) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus.x")
rownames(sig_OTUs_RFTM18)
#STILL SEQUENCE 



## SETTING OTU, sig_RF_TAX18, SAMP ####
OTU18 <- otu_table(otumat_matrix18, taxa_are_rows = FALSE)
sig_RF_TAX18 <- tax_table(sig_OTUs_RFTM18)
SAMP18 <- sample_data(metadata18_df)



## INSPECTING SAMPLE NAMES####
sample_names(SAMP18)
sample_names(OTU18)
sample_names(sig_RF_TAX18)

dim(SAMP18)
dim(OTU18)
dim(sig_RF_TAX18)


## EVENING OUT THE DATA ####
OTU18=transform_sample_counts(OTU18, function(x) 1E6 * x/sum(x))



## READING THROUGH PHYLOSEQ ####
physeq_sig18 = phyloseq(OTU18, sig_RF_TAX18, SAMP18) 
physeq_sig18


rftmdata18.ord <- ordinate(physeq_sig18, "NMDS", "bray")



view(SAMP18)
dim(SAMP18)
dim(OTU18)
dim(sig_RF_TAX18)



##NMDS - SIG OTUS - KINGDOM ####
plot_ordination(physeq_sig18, data_sig18.ord, type="taxa", color="Kingdom")+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#2E86C1"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS: RFTM OTUs - Kingdom",
             subtitle = "TAXA - NMDS Plot of Significant OTUs by RFTM",
             caption = "Data source: Oyster 16s 2018")

ggsave(filename = "NMDS_ROTU_Kingdom18.jpeg", plot=last_plot(), path ="Plot Diagrams 2018", width = 7, height = 5)  
 

##NMDS - SIG OTUS - PHYLUM ####
plot_ordination(physeq_sig18, data_sig18.ord, type="taxa", color="Phylum")+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#2E86C1"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS: RFTM OTUs - Phylum",
       subtitle = "TAXA - NMDS Plot of Significant OTUs by RFTM",
       caption = "Data source: Oyster 16s 2018")

ggsave(filename = "NMDS_ROTU_Phylum18.jpeg", plot=last_plot(), path ="Plot Diagrams 2018", width = 7, height = 5)  



##NMDS - SIG OTUS - CLASS ####
plot_ordination(physeq_sig18, data_sig18.ord, type="taxa", color="Class")+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#2E86C1"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS: RFTM OTUs - Class",
       subtitle = "TAXA - NMDS Plot of Significant OTUs by RFTM",
       caption = "Data source: Oyster 16s 2018")

ggsave(filename = "NMDS_ROTU_Class18.jpeg", plot=last_plot(), path ="Plot Diagrams 2018", width = 7, height = 5)  



##NMDS - SIG OTUS - ORDER ####
plot_ordination(physeq_sig18, data_sig18.ord, type="taxa", color="Order")+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#2E86C1"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS: RFTM OTUs - Order",
       subtitle = "TAXA - NMDS Plot of Significant OTUs by RFTM",
       caption = "Data source: Oyster 16s 2018")

ggsave(filename = "NMDS_ROTU_Order18.jpeg", plot=last_plot(), path ="Plot Diagrams 2018", width = 10, height = 5)   



##NMDS - SIG OTUS - FAMILY ####
plot_ordination(physeq_sig18, data_sig18.ord, type="taxa", color="Family")+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#2E86C1"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS: RFTM OTUs - Family",
       subtitle = "TAXA - NMDS Plot of Significant OTUs by RFTM",
       caption = "Data source: Oyster 16s 2018")

ggsave(filename = "NMDS_ROTU_Family18.jpeg", plot=last_plot(), path ="Plot Diagrams 2018", width = 10, height = 5)  


##NMDS - SIG OTUS - GENUS ####
plot_ordination(physeq_sig18, data_sig18.ord, type="taxa", color="Genus.x")+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#2E86C1"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS: RFTM OTUs - Genus",
       subtitle = "TAXA - NMDS Plot of Significant OTUs by RFTM",
       caption = "Data source: Oyster 16s 2018")

ggsave(filename = "NMDS_ROTU_Genus18.jpeg", plot=last_plot(), path ="Plot Diagrams 2018", width = 7, height = 5)  




#NMDS Plots for XXXX; type = "Samples"####

##NMDS - SIG OTUS - Species.x #### 
plot_ordination(physeq_sig18, data_sig18.ord, type="samples", color="Species.x")+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#2E86C1"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS: RFTM OTUs - Bivalve Species",
       subtitle = "SAMPLES - NMDS Plot of Significant OTUs by RFTM",
       caption = "Data source: Oyster 16s 2018", col="Bivalve Species")
  
ggsave(filename = "NMDS_ROTU_BiSpecies18.jpeg", plot=last_plot(), path ="Plot Diagrams 2018", width = 7, height = 5)  


##NMDS - SIG OTUS - Treatment2_18 #### 
plot_ordination(physeq_sig18, data_sig18.ord, type="samples", color="Treatment2_18")+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#2E86C1"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS: RFTM OTUs - Treatment Type",
       subtitle = "SAMPLES - NMDS Plot of Significant OTUs by RFTM",
       caption = "Data source: Oyster 16s 2018", col="Treatment Type")

ggsave(filename = "NMDS_ROTU_Treat2_18.jpeg", plot=last_plot(), path ="Plot Diagrams 2018", width = 7, height = 5)  


##NMDS - SIG OTUS - Treatment2_18 #### 
plot_ordination(physeq_sig18, data_sig18.ord, type="samples", color="Treatment2_18",)+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#2E86C1"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS: RFTM OTUs - Treatment Type",
       subtitle = "SAMPLES - NMDS Plot of Significant OTUs by RFTM",
       caption = "Data source: Oyster 16s 2018", col="Treatment Type")+
  facet_wrap(~ Species.x, )

ggsave(filename = "NMDS_ROTU_Treat2_Spec18.jpeg", plot=last_plot(), path ="Plot Diagrams 2018", width = 7, height = 5)  




#NMDS Plots for XXXX; type = "biplot"####
plot_ordination(physeq_sig18, data_sig18.ord, type="biplot", color = "Treatment2_18", shape ="Species.x" )+
  labs(col="Bivalve Species", type ="horsepower", shape = "SHAPE") 

 #Q: is biplot nessesary? It is hard to view because there is so much going on.
  





#Ledgend title manipulation #####
#labs(x="miles per gallon", y="displacement", size="horsepower", col="# of cylinders", shape="# of gears")

plot_bar(physeq_class18,"RFTM_score.x", fill="Species.x")+
  geom_col()+
  labs(title = "Barplot Site Frequency ",
       subtitle = "By RFTM Score",
       caption = "Data source: Oyster 16s 2017")+
  theme(legend.position="right", legend.text=element_text(size=10),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12),
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        panel.grid.major = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5))




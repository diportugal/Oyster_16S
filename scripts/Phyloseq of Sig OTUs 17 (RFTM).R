#Oyster 16S - Phyloseq Anaysis for 2017 RFTM Significant Data
#Author: Diana Portugal 
#Contact: dportugal8@gmail.com 
#Date 24-7-2021


#Installing packages####
library("phyloseq"); packageVersion("phyloseq")
data("GlobalPatterns")
library("ggplot2"); packageVersion("ggplot2")
library("plyr"); packageVersion("plyr")
library("data.table")
library("vegan")
library("dplyr")
library("ggpubr")
#install.packages("ggpubr")




#Isolating Significant OTUs for RFTM Variable (2017) ####

#sig_OTUs_RFTM17 <- RFTM_sig17 %>% 
  #select(Kingdom, Phylum, Class, Order, Family, Genus.x)

rftmpos = RFTM_sig17[RFTM_sig17$log2FoldChange<0,]
rftmasv = row.names(rftmpos)

phylo17_rftm = prune_taxa(taxa_sums(physeq_class) > 0,physeq_class ) 

rftmdata17.ord <- ordinate(phylo17_rftm, "NMDS", "bray")




#ALPHA DIVERSITY ####
#The plots don't change with the different phyloseq objects

plot_richness(physeq_class, measures=c("Chao1", "Shannon"), x="RFTM_score.x", col="Site.x")

plot_richness(physeq_sigPea, measures=c("Chao1", "Shannon"), x="peacrabs.x", col="Site.x")

plot_richness(RFTM_sig17, measures=c("Chao1", "Shannon"), x="peacrabs.x", col="Site.x")



##Site ####
comparison_SigRtfm_alph <- list(c("NW", "OY"), c("OY", "SW"), c("NW", "SW"))
symnum.args17_w = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

plot_richness(physeq_class, x="Site.x", measures=c("Chao1", "Shannon"), color = "Site.x")+
  geom_boxplot(alpha=0.6)+
  labs(title = "Alpha Diversity",
       subtitle = "By Site",
       caption = "Data source: Oyster 16s 2017")+
  theme(legend.position="right", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12),
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        panel.grid.major = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5))+
  stat_compare_means(method = "wilcox.test", comparisons = comparison_SigRtfm_alph, label = "p.signif") 
ggsave(filename = "AD_By_Site_17.jpeg", plot=last_plot(), path ="Plot Diagrams/Alpha Diversity/", width = 7, height = 5)  


plot_richness(physeq_class, x="Site.x", measures=c("Chao1", "Shannon"), color = "Site.x")+
  geom_boxplot(alpha=0.6)+
  labs(title = "Alpha Diversity",
       subtitle = "By Site & RFTM Score",
       caption = "Data source: Oyster 16s 2017")+
  theme(legend.position="right", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12),
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        panel.grid.major = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5))+
  stat_compare_means(method = "wilcox.test", comparisons = comparison_SigRtfm_alph, label = "p.signif")+
  facet_wrap(~RFTM_score.x, 2)
ggsave(filename = "AD_BySite_ByRFTM_17.jpeg", plot=last_plot(), path ="Plot Diagrams/Alpha Diversity/", width = 7, height = 5)  


plot_richness(physeq_class, x="Site.x", measures=c("Chao1", "Shannon"), color = "Site.x")+
  geom_boxplot(alpha=0.6)+
  labs(title = "Alpha Diversity",
       subtitle = "By Site & Pea Crab",
       caption = "Data source: Oyster 16s 2017")+
  theme(legend.position="right", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12),
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        panel.grid.major = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5))+
  stat_compare_means(method = "wilcox.test", comparisons = comparison_SigRtfm_alph, label = "p.signif")+
  facet_wrap(~peacrabs.x, 2)
ggsave(filename = "AD_BySite_ByPea.jpeg", plot=last_plot(), path ="Plot Diagrams/Alpha Diversity/", width = 7, height = 5)  





##Treatment2 ####
comparison_treat_alph <- list(c("HIGH_MONO", "HIGH_POLY"), c("HIGH_MONO", "LOW_MONO"), c("HIGH_MONO", "LOW_POLY"), c("HIGH_POLY", "LOW_MONO"), c("HIGH_POLY", "LOW_POLY"),c("LOW_MONO", "LOW_POLY"))
symnum.args17_w = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("*****","****", "***", "**", "*", "ns"))

plot_richness(physeq_class, x="Treatment2", measures=c("Chao1", "Shannon"), color = "Treatment2")+
  geom_boxplot(alpha=0.6)+
  labs(title = "Alpha Diversity",
       subtitle = "By Treatment",
       caption = "Data source: Oyster 16s 2017")+
  theme(legend.position="right", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12),
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        panel.grid.major = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5))+
  stat_compare_means(method = "wilcox.test", comparisons = comparison_treat_alph, label = "p.signif")
ggsave(filename = "AD_ByTreat_17.jpeg", plot=last_plot(), path ="Plot Diagrams/Alpha Diversity/", width = 7, height = 5)  


plot_richness(physeq_class, x="Treatment2", measures=c("Chao1", "Shannon"), color = "Treatment2")+
  geom_boxplot(alpha=0.6)+
  labs(title = "Alpha Diversity",
       subtitle = "By Treatment & RFTM Score",
       caption = "Data source: Oyster 16s 2017")+
  theme(legend.position="right", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12),
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        panel.grid.major = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5))+
  stat_compare_means(method = "wilcox.test", comparisons = comparison_treat_alph, label = "p.signif")+
  facet_wrap(~RFTM_score.x, 2) 
ggsave(filename = "AD_ByTreat_ByRFTM_17.jpeg", plot=last_plot(), path ="Plot Diagrams/Alpha Diversity/", width = 7, height = 10)  


plot_richness(physeq_class, x="Treatment2", measures=c("Chao1", "Shannon"), color = "Treatment2")+
  geom_boxplot(alpha=0.6)+
  labs(title = "Alpha Diversity",
       subtitle = "By Treatment & Pea Crab",
       caption = "Data source: Oyster 16s 2017")+
  theme(legend.position="right", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12),
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        panel.grid.major = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5))+
  stat_compare_means(method = "wilcox.test", comparisons = comparison_treat_alph, label = "p.signif")+
  facet_wrap(~peacrabs.x,)
ggsave(filename = "AD_ByTreat_Pea_17.jpeg", plot=last_plot(), path ="Plot Diagrams/Alpha Diversity/", width = 7, height = 7)  


plot_richness(physeq_class, x="RFTM_score.x", measures=c("Chao1", "Shannon"), color = "RFTM_score.x")+
  labs(title = "Alpha Diversity",
       subtitle = "By Treatment & Pea Crab",
       caption = "Data source: Oyster 16s 2017")+
  theme(legend.position="right", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12),
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        panel.grid.major = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5))
  



# BAR PLOT ####
#comment 

## Site Frequency by RFTM Score ####

plot_bar(physeq_class,"RFTM_score.x", fill="Null")+
  geom_col()+
  scale_colour_manual(values=c("#FA7169","#2E86C1","#8FC172","#AF7AC5"))+
  labs(title = "Barplot Site Frequency ",
       subtitle = "By RFTM Score",
       caption = "Data source: Oyster 16s 2017")+
  theme(legend.position="right", legend.text=element_text(size=10),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12),
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        panel.grid.major = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5))+
  facet_wrap(~Site.x,)
ggsave(filename = "BP_RFTM_17.jpeg", plot=last_plot(), path ="Plot Diagrams/Bar Plot/", width = 7, height = 7)  



plot_bar(physeq_class,"RFTM_score.x", fill="Treatment2")+
  geom_col()+
  scale_colour_manual(values=c("#FA7169","#2E86C1","#8FC172","#AF7AC5" ))+
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
ggsave(filename = "BP_RFTM_17.jpeg", plot=last_plot(), path ="Plot Diagrams/Bar Plot/", width = 7, height = 7)  






plot_bar(physeq_class,"peacrabs.x", fill="Site.x")+
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







#Phyloseq of Sig PeaCrab ####

#Isolating significant OTUs for RFTM Score
phylo_peacrab <- PEA_sig17 %>% 
  select(Kingdom, Phylum, Class, Order, Family, Genus.x)


metadata17 <- read.csv("Oyster_data_raw/meta17cleaned")
asvtable17 <- fread("Oyster_data_raw/asvtable_de17.csv")
run23 <- read.csv("Oyster_data_raw/Run123_taxa_complete.csv")


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

## CONVERTING TO MATRICIES
metadata17_df <- as.data.frame(metadata17, rownames("X.1"))
rownames(metadata17_df)
#STILL UNIQUE ID

otumat_matrix <- as.matrix(asvtable17, rownames=rownames(asvtable17))
rownames(otumat_matrix)
#STILL UNIQUE ID 

#Using significant otu data from DESeq2 analysis of RFTM score
sig_OTUs_peacrab <- as.matrix(phylo_peacrab)
colnames(phylo_peacrab) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus.x")
rownames(sig_OTUs_peacrab)
#STILL SEQUENCE 



## SETTING OTU, sig_RF_TAX18, SAMP
OTU <- otu_table(otumat_matrix, taxa_are_rows = FALSE)
sig_PEA_TAX <- tax_table(sig_OTUs_peacrab)
SAMP <- sample_data(metadata17_df)



## INSPECTING SAMPLE NAMES
sample_names(SAMP)
sample_names(OTU)
sample_names(sig_PEA_TAX)

dim(SAMP)
dim(OTU)
dim(sig_PEA_TAX)


## EVENING OUT THE DATA
OTU=transform_sample_counts(OTU, function(x) 1E6 * x/sum(x))

## READING THROUGH PHYLOSEQ
physeq_sigPea = phyloseq(OTU, sig_PEA_TAX, SAMP) 
physeq_sigPea

peacrabdata.ord <- ordinate(physeq_sigPea, "NMDS", "bray")


































































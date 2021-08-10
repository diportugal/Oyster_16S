#NMDS IN DEPTH 
#2021-06-29
#Author:Diana Portugal
#Contact: dportugal8@gmail.com 



library("phyloseq"); packageVersion("phyloseq")
data("GlobalPatterns")
library("ggplot2"); packageVersion("ggplot2")
library("plyr"); packageVersion("plyr")
theme_set(theme_bw())


#ORIGINAL DATA ####
physeq_class


#FILTERING FOR THE TOP 5 TAXA####

## Kingdom ####
#Filtering for top 5 Kingdom (NOT DONE, not necessary)

## Phylum ####
#filtering for top 5 phyla 
phylum.sum = tapply(taxa_sums(physeq_class), tax_table(physeq_class)[, "Phylum"], sum, na.rm=TRUE)
top5phyla_nmds = names(sort(phylum.sum, TRUE))[1:5]
physeq_class_5phyla = prune_taxa((tax_table(physeq_class)[, "Phylum"] %in% top5phyla_nmds), physeq_class)


## Class ####
#Filtering for top 5 Class 
class.sum = tapply(taxa_sums(physeq_class), tax_table(physeq_class)[, "Class"], sum, na.rm=TRUE)
top5class_nmds = names(sort(class.sum, TRUE))[1:5]
physeq_class_5class = prune_taxa((tax_table(physeq_class)[, "Class"] %in% top5class_nmds), physeq_class)


## Order ####
#Filtering for top 5 Order 
order.sum = tapply(taxa_sums(physeq_class), tax_table(physeq_class)[, "Order"], sum, na.rm=TRUE)
top5order_nmds = names(sort(order.sum, TRUE))[1:5]
physeq_class_5order = prune_taxa((tax_table(physeq_class)[, "Order"] %in% top5order_nmds), physeq_class)


## Family ####
#Filtering for top 5 Family
family.sum = tapply(taxa_sums(physeq_class), tax_table(physeq_class)[, "Family"], sum, na.rm=TRUE)
top5family_nmds = names(sort(family.sum, TRUE))[1:5]
physeq_class_5family = prune_taxa((tax_table(physeq_class)[, "Family"] %in% top5family_nmds), physeq_class)


## Genus ####
#Filtering for top 5 genus 
genus.sum = tapply(taxa_sums(physeq_class), tax_table(physeq_class)[, "Genus.x"], sum, na.rm=TRUE)
top5genus_nmds = names(sort(genus.sum, TRUE))[1:5]
physeq_class_5genus = prune_taxa((tax_table(physeq_class)[, "Genus.x"] %in% top5genus_nmds), physeq_class)



#REMOVING NA(s) FROM EACH LEVEL ####

## No Kingdom ####

## Phylum ####
physeq_class_5phyla = prune_taxa((tax_table(physeq_class)[, "Phylum"] %in% top5phyla_nmds), physeq_class)
physeq_phylum5_no_na <- subset_samples(physeq_class_5phyla, "Phylum" != "NA", )
data17.ord_phylum5_no_na <- ordinate(physeq_phylum5_no_na, "NMDS", "bray")

## Class ####
physeq_class_5class = prune_taxa((tax_table(physeq_class)[, "Class"] %in% top5class_nmds), physeq_class)
physeq_class5_no_na <- subset_samples(physeq_class_5class, "Class" != "NA", )
data17.ord_class5_no_na <- ordinate(physeq_class5_no_na, "NMDS", "bray")

## Order ####
physeq_class_5order = prune_taxa((tax_table(physeq_class)[, "Order"] %in% top5order_nmds), physeq_class)
physeq_order5_no_na <- subset_samples(physeq_class_5order, "Order" != "NA", )
data17.ord_order5_no_na <- ordinate(physeq_order5_no_na, "NMDS", "bray")

## Family ####
physeq_class_5family = prune_taxa((tax_table(physeq_class)[, "Family"] %in% top5family_nmds), physeq_class)
physeq_family5_no_na <- subset_samples(physeq_class_5family, "Family" != "NA", )
data17.ord_family5_no_na <- ordinate(physeq_family5_no_na, "NMDS", "bray")

## Genus.x ####
physeq_class_5genus = prune_taxa((tax_table(physeq_class)[, "Genus.x"] %in% top5genus_nmds), physeq_class)
physeq_genus5_no_na <- subset_samples(physeq_class_5genus, "Genus.x" != "NA", )
data17.ord_genus5_no_na <- ordinate(physeq_genus5_no_na, "NMDS", "bray")



#PLOT MANIPULATION EXAMPLES ####
#STAGE 1 (Basic form)
data17.ord <- ordinate(physeq_class, "NMDS", "bray")
nmds_plot1= plot_ordination(physeq_class_5phyla, data17.ord, type="taxa", color="Phylum", title="NMDS of Top 5 Phyla")
print(nmds_plot1)

#no na 
plot_ordination(physeq_phylum5_no_na, data17.ord_phylum5_no_na, type="taxa", color="Phylum", title="NMDS of Top 5 Phyla")


#STAGE 2 (Color Manipulation)
nmds_plot2= plot_ordination(physeq_class_5phyla, data17.ord, type="taxa", color="Phylum", title="NMDS of Top 5 Phyla")+
  scale_colour_manual(values=c("Cyanobacteria" = "#FA7169", "Firmicutes" = "#33BEFF", "Proteobacteria" = "#8FC172", "Spirochaetota"= "#FAC069", "Verrucomicrobiota"="#BA96D9", "NA"= NULL )) +
  coord_equal() +
  theme_bw()
print(nmds_plot2)
#This prints out the difference in the ordination table of the taxa calculated based on the top 5 phyla found 
#Need to null out the NA value that is GREY


#STAGE 3 (Legend Manipulation)
nmds_plot3= plot_ordination(physeq_class_5phyla, data17.ord, type="taxa", color="Phylum", title="NMDS of Top 5 Phyla")+
  scale_colour_manual(values=c("Cyanobacteria" = "#FA7169", "Firmicutes" = "#33BEFF", "Proteobacteria" = "#8FC172", "Spirochaetota"= "#FAC069", "Verrucomicrobiota"="#BA96D9", "NA"= NULL ))+
  theme(legend.position="left", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10))
print(nmds_plot3)
#We can put the legend on the "top", "bottom", "right", or "left" 



#STAGE 4 (Title Manipulation)
#NMDS Phyla
nmds_plot4= plot_ordination(physeq_class_5phyla, data17.ord, type="taxa", color="Phylum")+
  scale_colour_manual(values=c("Cyanobacteria" = "#FA7169", "Firmicutes" = "#33BEFF", "Proteobacteria" = "#8FC172", "Spirochaetota"= "#FAC069", "Verrucomicrobiota"="#BA96D9", "NA"= NULL ))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS of Top 5 Phyla",
       subtitle = "Plot the top 5 phyla in Oyster 16s data",
       caption = "Data source: Oyster 16s 2017")
print(nmds_plot4)
#We can bold the title, adjust to the left, center, and right, manipulate the size, and color
#plot.title=element_text(hjust = )
#0 = left, 0.5 = center, 1 = right 
#We can add a subtitle and data citation for figures






#ANALYSIS BY TAXA ####
#MAKING NMDS (TAXA) PLOTS FOR THE TOP RESULTS ON EACH TAXA 
#Only considering taxa - without NA values

#NMDS of Kingdom Plots 
plot_ordination(physeq_class, data17.ord, type="taxa", color="Kingdom")+
  scale_colour_manual(values=c("Archaea" = "#FA7169", "Bacteria" = "#2E86C1", "Eukaryota" = "#8FC172", "NA"= "#AF7AC5" ))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS Kingdom Plots",
       subtitle = "Plotting Kingdom Diversity (Top 5)",
       caption = "Data source: Oyster 16s 2017")
ggsave(filename = "NMDS_Top5King_17.jpeg", plot=last_plot(), path ="Plot Diagrams/NMDS Original /", width = 7, height = 5)  


#NMDS of Phylum Plots 
plot_ordination(physeq_phylum5_no_na, data17.ord_phylum5_no_na, type="taxa", color="Phylum")+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS Phylum Plots",
       subtitle = "Plotting Phyla Diversity (Top 5) without NAs",
       caption = "Data source: Oyster 16s 2017")
ggsave(filename = "NMDS_Top5Phyla_17.jpeg", plot=last_plot(), path ="Plot Diagrams/NMDS Original /", width = 7, height = 5)  


#NMDS of Class 
plot_ordination(physeq_class5_no_na, data17.ord_class5_no_na, type="taxa", color="Class")+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS Class Plots",
       subtitle = "Plotting Class Diversity (Top 5) without NAs",
       caption = "Data source: Oyster 16s 2017")
ggsave(filename = "NMDS_Top5Class_17.jpeg", plot=last_plot(), path ="Plot Diagrams/NMDS Original /", width = 7, height = 5)  



#NMDS of Order 
plot_ordination(physeq_order5_no_na, data17.ord_order5_no_na, type="taxa", color="Order")+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS Order Plots",
       subtitle = "Plotting Order Diversity (Top 5) without NAs",
       caption = "Data source: Oyster 16s 2017")
ggsave(filename = "NMDS_Top5Order_17.jpeg", plot=last_plot(), path ="Plot Diagrams/NMDS Original /", width = 7, height = 5)  
  

#NMDS of Family 
plot_ordination(physeq_family5_no_na, data17.ord_family5_no_na, type="taxa", color="Family")+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS Family Plots",
       subtitle = "Plotting Family Diversity (Top 5) without NAs",
       caption = "Data source: Oyster 16s 2017")
ggsave(filename = "NMDS_Top5Family_17.jpeg", plot=last_plot(), path ="Plot Diagrams/NMDS Original /", width = 7, height = 5)  


#NMDS of Genus 
plot_ordination(physeq_genus5_no_na, data17.ord_genus5_no_na, type="taxa", color="Genus.x")+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS Genus Plots",
       subtitle = "Plotting Genus Diversity (Top 5) without NAs",
       caption = "Data source: Oyster 16s 2017")
ggsave(filename = "NMDS_Top5Genus_17.jpeg", plot=last_plot(), path ="Plot Diagrams/NMDS Original /", width = 7, height = 5)  



#ANALYSIS BY SAMPLES ####

#This is showing that the oysters that have a peacrab present a less similar than the ones that have no peacrab presence 
plot_ordination(physeq_class, data17.ord, type="samples", color="peacrabs.x", shape=NULL) + 
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS Peacrab Sample Plot",
       subtitle = "Plotting Peacrab Sample Diversity",
       caption = "Data source: Oyster 16s 2017")



plot_ordination(physeq_class, data17.ord, type="samples", color="RFTM_score.x", shape=NULL)+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS RFTM Score Sample Plot",
       subtitle = "Plotting RFTM Score Diversity",
       caption = "Data source: Oyster 16s 2017")+
  facet_wrap(~Site.x, )


plot_ordination(physeq_class, data17.ord, type="samples", color="RFTM_score.x", shape="Site.x")+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS RFTM/Site/Peacrab",
       subtitle = "Community Diversity",
       caption = "Data source: Oyster 16s 2017")+
  facet_wrap(~peacrabs.x, )




#The next two plots are the same but are facet wrapped by different variables:

#Treatment2 by RFTM Score separated by RFTM Scores
plot_ordination(physeq_class, data17.ord, type="samples", color="RFTM_score.x", shape="Treatment2")+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS Treatment/RFTM Score Sample Plot",
       subtitle = "Plotting Treatment/RFTM Score Diversity",
       caption = "Data source: Oyster 16s 2017")+
  facet_wrap(~RFTM_score.x, 2)


#Treatment2 by RFTM Score separated by Treatment2
plot_ordination(physeq_class, data17.ord, type="samples", color="RFTM_score.x", shape="Treatment2")+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS RFTM Score/Treatment Sample Plot",
       subtitle = "Plotting RFTM Score/Treatment Diversity",
       caption = "Data source: Oyster 16s 2017")



plot_ordination(physeq_class, data17.ord, type="samples", color="RFTM_score.x", shape=NULL)+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS RFTM Score/Treatment Sample Plot",
       subtitle = "Plotting RFTM Score/Treatment Diversity",
       caption = "Data source: Oyster 16s 2017")+
  facet_wrap(~peacrabs.x, )



#creating labels for peacrabs ####
plot_ordination(physeq_class, data17.ord, type="samples", color="RFTM_score.x", shape="Site.x", label="X")+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS RFTM/Site/Peacrab",
       subtitle = "Community Diversity",
       caption = "Data source: Oyster 16s 2017")+
  facet_wrap(~peacrabs.x, )










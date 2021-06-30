#NMDS IN DEPTH 
#2021-06-29
#Author:Diana Portugal
#Contact: dportugal8@gmail.com 



library("phyloseq"); packageVersion("phyloseq")
data("GlobalPatterns")
library("ggplot2"); packageVersion("ggplot2")
library("plyr"); packageVersion("plyr")
theme_set(theme_bw())


#Filtering for top 5 Kingdom (NOT DONE, not really necessary)
genus.sum = tapply(taxa_sums(physeq_class), tax_table(physeq_class)[, "Genus"], sum, na.rm=TRUE)
top5genus_nmds = names(sort(genus.sum, TRUE))[1:5]
physeq_class_5genus = prune_taxa((tax_table(physeq_class)[, "Genus"] %in% top5genus_nmds), physeq_class)


#filtering for top 5 phyla 
phylum.sum = tapply(taxa_sums(physeq_class), tax_table(physeq_class)[, "Phylum"], sum, na.rm=TRUE)
top5phyla_nmds = names(sort(phylum.sum, TRUE))[1:5]
physeq_class_5phyla = prune_taxa((tax_table(physeq_class)[, "Phylum"] %in% top5phyla_nmds), physeq_class)


#Filtering for top 5 Class (NOT DONE)
genus.sum = tapply(taxa_sums(physeq_class), tax_table(physeq_class)[, "Genus"], sum, na.rm=TRUE)
top5genus_nmds = names(sort(genus.sum, TRUE))[1:5]
physeq_class_5genus = prune_taxa((tax_table(physeq_class)[, "Genus"] %in% top5genus_nmds), physeq_class)


#Filtering for top 5 Order (NOT DONE)
genus.sum = tapply(taxa_sums(physeq_class), tax_table(physeq_class)[, "Genus"], sum, na.rm=TRUE)
top5genus_nmds = names(sort(genus.sum, TRUE))[1:5]
physeq_class_5genus = prune_taxa((tax_table(physeq_class)[, "Genus"] %in% top5genus_nmds), physeq_class)


#Filtering for top 5 Family (NOT DONE)
genus.sum = tapply(taxa_sums(physeq_class), tax_table(physeq_class)[, "Genus"], sum, na.rm=TRUE)
top5genus_nmds = names(sort(genus.sum, TRUE))[1:5]
physeq_class_5genus = prune_taxa((tax_table(physeq_class)[, "Genus"] %in% top5genus_nmds), physeq_class)


#Filtering for top 5 genus 
genus.sum = tapply(taxa_sums(physeq_class), tax_table(physeq_class)[, "Genus"], sum, na.rm=TRUE)
top5genus_nmds = names(sort(genus.sum, TRUE))[1:5]
physeq_class_5genus = prune_taxa((tax_table(physeq_class)[, "Genus"] %in% top5genus_nmds), physeq_class)

#WE COULD FILTER OUT THE TOP 5 RESULTS FOR EACH TAXA (IF WE WANTED TO....)





#STAGE 1 (Basic form)
data17.ord <- ordinate(physeq_class, "NMDS", "bray")
nmds_plot1= plot_ordination(physeq_class_5phyla, data17.ord, type="taxa", color="Phylum", title="NMDS of Top 5 Phyla")
print(nmds_plot1)


#STAGE 2 (Color Manipulation)
data17.ord <- ordinate(physeq_class, "NMDS", "bray")
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
nmds_plot4= plot_ordination(physeq_class_5phyla, data17.ord, type="taxa", color="Phylum")+
  scale_colour_manual(values=c("Cyanobacteria" = "#FA7169", "Firmicutes" = "#33BEFF", "Proteobacteria" = "#8FC172", "Spirochaetota"= "#FAC069", "Verrucomicrobiota"="#BA96D9", "NA"= NULL ))+
  theme(legend.position="left", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 20, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS of Top 5 Phyla",
       subtitle = "Plot the top 5 phyla in Oyster 16s data",
       caption = "Data source: Oyster 16s 2017")
print(nmds_plot4)
#We can bold the title, adjust to the left, center, and right, manipulate the size, and color
#plot.title=element_text(hjust = )
#0 = left, 0.5 = center, 1 = right 
#We can add a subtitle and data citation for figures


nmds_plot4 + facet_wrap(~Phylum, 3)
#Separating out the 5 different Phyla + NA
#Stack of 3


#We can pick variables from only the metadata17_df table that has the information of each oyster sample
nmds_samp1 = plot_ordination(physeq_class_5phyla, data17.ord, type="samples", color="Treatment2", shape="Color") 
nmds_samp1 + geom_polygon(aes(fill=Site.x)) + geom_point(size=5) + ggtitle(":)")































#Oyster 16S NMDS Anaysis for 2018 Data
#Author: Diana Portugal 
#Contact: dportugal8@gmail.com 



library("phyloseq"); packageVersion("phyloseq")
data("GlobalPatterns")
library("ggplot2"); packageVersion("ggplot2")
library("plyr"); packageVersion("plyr")
theme_set(theme_bw())


#ORIGINAL DATA ####
physeq_class18


data18.ord <- ordinate(physeq_class18, "NMDS", "bray")




#FILTERING FOR THE TOP 5 TAXA ####

#Filtering for top 5 Kingdom  (NOT DONE, not really necessary)
kingdom.sum18 = tapply(taxa_sums(physeq_class18), tax_table(physeq_class18)[, "Kingdom"], sum, na.rm=TRUE)
top5kingdom_nmds18 = names(sort(kingdom.sum18, TRUE))[1:5]
physeq_class_5kingdom18 = prune_taxa((tax_table(physeq_class18)[, "Kingdom"] %in% top5kingdom_nmds18), physeq_class18)


#filtering for top 5 phyla  
phylum.sum18 = tapply(taxa_sums(physeq_class18), tax_table(physeq_class18)[, "Phylum"], sum, na.rm=TRUE)
top5phyla_nmds18 = names(sort(phylum.sum18, TRUE))[1:5]
physeq_class_5phyla18 = prune_taxa((tax_table(physeq_class18)[, "Phylum"] %in% top5phyla_nmds18), physeq_class18)


#Filtering for top 5 Class 
class.sum18 = tapply(taxa_sums(physeq_class18), tax_table(physeq_class18)[, "Class"], sum, na.rm=TRUE)
top5class_nmds18 = names(sort(class.sum18, TRUE))[1:5]
physeq_class_5class18 = prune_taxa((tax_table(physeq_class18)[, "Class"] %in% top5class_nmds18), physeq_class18)


#Filtering for top 5 Order 
order.sum18 = tapply(taxa_sums(physeq_class18), tax_table(physeq_class18)[, "Order"], sum, na.rm=TRUE)
top5order_nmds18 = names(sort(order.sum18, TRUE))[1:5]
physeq_class_5order18 = prune_taxa((tax_table(physeq_class18)[, "Order"] %in% top5order_nmds18), physeq_class18)


#Filtering for top 5 Family 
family.sum18 = tapply(taxa_sums(physeq_class18), tax_table(physeq_class18)[, "Family"], sum, na.rm=TRUE)
top5family_nmds18 = names(sort(family.sum18, TRUE))[1:5]
physeq_class_5family18 = prune_taxa((tax_table(physeq_class18)[, "Family"] %in% top5family_nmds18), physeq_class18)


#Filtering for top 5 genus
genus.sum18 = tapply(taxa_sums(physeq_class18), tax_table(physeq_class18)[, "Genus.x"], sum, na.rm=TRUE)
top5genus_nmds18 = names(sort(genus.sum18, TRUE))[1:5]
physeq_class_5genus18 = prune_taxa((tax_table(physeq_class18)[, "Genus.x"] %in% top5genus_nmds18), physeq_class18)




#Evaluation of Top 5 Taxa ####


#NMDS Plot evaluating 2018 data Kingdom Taxa
plot_ordination(physeq_class18, data18.ord, type="taxa", color="Kingdom")+
  scale_colour_manual(values=c("Archaea" = "#FA7169", "Bacteria" = "#2E86C1", "Eukaryota" = "#8FC172", "NA"= "#AF7AC5" ))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS Kingdom Plots",
       subtitle = "Plotting Top 5 Kingdoms Diversity",
       caption = "Data source: Oyster 16s 2018")

ggsave(filename = "NMDS_Tax_King18.jpeg", plot=last_plot(), path ="Plot Diagrams 2018/Original_2018/", width = 7, height = 5)  



#NMDS Plot evaluating 2018 data Phylum Taxa
plot_ordination(physeq_class_5phyla18, data18.ord, type="taxa", color="Phylum")+
  scale_colour_manual(values=c("Bacteroidota" = "#FA7169", "Cyanobacteria" = "#2E86C1", "Proteobacteria" = "#8FC172", "Spirochaetota"= "#AF7AC5", "Verrucomicrobiota"="#3FD8FF"))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS Phylum Plots",
       subtitle = "Plotting Top 5 Phylum Diversity",
       caption = "Data source: Oyster 16s 2018")

ggsave(filename = "NMDS_Tax_TopPhy18.jpeg", plot=last_plot(), path ="Plot Diagrams 2018/Original_2018/", width = 7, height = 5)  




#NMDS Plot evaluating 2018 data Class Taxa
plot_ordination(physeq_class_5class18, data18.ord, type="taxa", color="Class")+
  scale_colour_manual(values=c("Alphaproteobacteria" = "#FA7169", "Bacteroidia" = "#2E86C1", "Chlamydiae" = "#8FC172", "Gammaproteobacteria"= "#AF7AC5"))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS Class Plots",
       subtitle = "Plotting Top 5 Class Diversity",
       caption = "Data source: Oyster 16s 2018")

ggsave(filename = "NMDS_Tax_TopClass18.jpeg", plot=last_plot(), path ="Plot Diagrams 2018/Original_2018/", width = 7, height = 5)  


#NMDS Plot evaluating 2018 data Order Taxa (Not Done)
plot_ordination(physeq_class_5order18, data18.ord, type="taxa", color="Order")+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS Order Plots",
       subtitle = "Plotting Top 5 Order Diversity",
       caption = "Data source: Oyster 16s 2018")

ggsave(filename = "NMDS_Tax_TopOrder18.jpeg", plot=last_plot(), path ="Plot Diagrams 2018/Original_2018/", width = 7, height = 5)  



#NMDS Plot evaluating 2018 data Family Taxa 
plot_ordination(physeq_class_5family18, data18.ord, type="taxa", color="Family")+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS Family Plots",
       subtitle = "Plotting Top 5 Family Diversity",
       caption = "Data source: Oyster 16s 2018")

ggsave(filename = "NMDS_Tax_TopFam18.jpeg", plot=last_plot(), path ="Plot Diagrams 2018/Original_2018/", width = 7, height = 5)  



#NMDS Plot evaluating 2018 data Genus Taxa
plot_ordination(physeq_class_5genus18, data18.ord, type="taxa", color="Genus.x")+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS Genus Plots",
       subtitle = "Plotting Top 5 Genus Diversity",
       caption = "Data source: Oyster 16s 2018")

ggsave(filename = "NMDS_Tax_TopGen18.jpeg", plot=last_plot(), path ="Plot Diagrams 2018/Original_2018/", width = 7, height = 5)  




##Evaluation of Samples ####

#Show the similarity between species 
plot_ordination(physeq_class18, data18.ord, type="samples", color="Species.x")+
  scale_colour_manual(values=c("AM" = "#FA7169", "CV" = "#AF7AC5", "IR" = "#8FC172", "LP"= "#2E86C1"))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#2E86C1"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS Species Differences",
       caption = "Data source: Oyster 16s 2018", col="Bivalve Species")

ggsave(filename = "NMDS_BiSpecies18.jpeg", plot=last_plot(), path ="Plot Diagrams 2018/Original_2018/", width = 7, height = 5)  



plot_ordination(physeq_class18, data18.ord, type="samples", color="Treatment2_18")+
 theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#2E86C1"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS Treatment Differences",
       caption = "Data source: Oyster 16s 2018", col="Treatment Type")

ggsave(filename = "NMDS_Treat2_18.jpeg", plot=last_plot(), path ="Plot Diagrams 2018/Original_2018/", width = 7, height = 5)  




plot_ordination(physeq_class18, data18.ord, type="samples", color="Treatment2_18")+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#2E86C1"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS Treatment Differences",
       caption = "Data source: Oyster 16s 2018", col="Treatment Type")+
  facet_wrap(~Species.x, )

ggsave(filename = "NMDS_Treat2_Spec18.jpeg", plot=last_plot(), path ="Plot Diagrams 2018/Original_2018/", width = 7, height = 5)  








#RFTM score has to be a factor here 
RFTM.num18 <- sample_data(physeq_class18)$RFTM.asnum18 = sample_data(physeq_class18)$RFTM.asnum18 = as.numeric(sample_data(physeq_class18)$RFTM_score.x)

plot_ordination(physeq_class18, data18.ord, type="samples", color = "RFTM_score.x")+
  scale_color_manual(values=c("#6495ED", "#FFAC1E", "#FF0000", "#02A026", "#525151", "#FF00FF", "#9700FF"))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#2E86C1"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "Community Differences by RFTM Score",
       caption = "Data source: Oyster 16s 2018")



plot_ordination(physeq_class18, data18.ord, type="samples", color = "RFTM_score.x")+
  scale_color_manual(values=c("#6495ED", "#FFAC1E", "#FF0000", "#02A026", "#525151", "#FF00FF", "#9700FF"))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#2E86C1"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "RFTM Score by Species",
       caption = "Data source: Oyster 16s 2018")+
  facet_wrap(~Species.x, )

 
#Busy Busy plot, not a lot of information derived from here - Don't do this quite yet
plot_ordination(physeq_class18, data18.ord, type="samples", color = "Weight")+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#2E86C1"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "yyy",
       caption = "Data source: Oyster 16s 2018")


#^^^^
#Question: How to combine the two data sets so that it doesnt create a Taxa and Samples shape and color ####




#Evaluation of Taxa/Samples ####

plot_ordination(physeq_class18, data18.ord, type="samples", color="delta_weight18")


plot_ordination(physeq_class18, data18.ord, type="biplot", color = "Species.x")




plot_ordination(physeq_class18, data18.ord, type="samples", color = "RFTM_score.x")+
  scale_color_manual(values=c("#6495ED", "#FFAC1E", "#FF0000", "#02A026", "#525151", "#FF00FF", "#9700FF"))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#2E86C1"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "Community Differences by RFTM Score",
       caption = "Data source: Oyster 16s 2018")+
  facet_wrap(~peacrab.x, )




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

##Kingdom #### 
#(NOT DONE, not really necessary)
kingdom.sum18 = tapply(taxa_sums(physeq_class18), tax_table(physeq_class18)[, "Kingdom"], sum, na.rm=TRUE)
top5kingdom_nmds18 = names(sort(kingdom.sum18, TRUE))[1:5]
physeq_class_5kingdom18 = prune_taxa((tax_table(physeq_class18)[, "Kingdom"] %in% top5kingdom_nmds18), physeq_class18)


##Phylum #### 
phylum.sum18 = tapply(taxa_sums(physeq_class18), tax_table(physeq_class18)[, "Phylum"], sum, na.rm=TRUE)
top5phyla_nmds18 = names(sort(phylum.sum18, TRUE))[1:5]
physeq_class_5phyla18 = prune_taxa((tax_table(physeq_class18)[, "Phylum"] %in% top5phyla_nmds18), physeq_class18)


##Class ####
class.sum18 = tapply(taxa_sums(physeq_class18), tax_table(physeq_class18)[, "Class"], sum, na.rm=TRUE)
top5class_nmds18 = names(sort(class.sum18, TRUE))[1:5]
physeq_class_5class18 = prune_taxa((tax_table(physeq_class18)[, "Class"] %in% top5class_nmds18), physeq_class18)


##Order #### 
order.sum18 = tapply(taxa_sums(physeq_class18), tax_table(physeq_class18)[, "Order"], sum, na.rm=TRUE)
top5order_nmds18 = names(sort(order.sum18, TRUE))[1:5]
physeq_class_5order18 = prune_taxa((tax_table(physeq_class18)[, "Order"] %in% top5order_nmds18), physeq_class18)


##Family ####
family.sum18 = tapply(taxa_sums(physeq_class18), tax_table(physeq_class18)[, "Family"], sum, na.rm=TRUE)
top5family_nmds18 = names(sort(family.sum18, TRUE))[1:5]
physeq_class_5family18 = prune_taxa((tax_table(physeq_class18)[, "Family"] %in% top5family_nmds18), physeq_class18)


##Genus ####
genus.sum18 = tapply(taxa_sums(physeq_class18), tax_table(physeq_class18)[, "Genus.x"], sum, na.rm=TRUE)
top5genus_nmds18 = names(sort(genus.sum18, TRUE))[1:5]
physeq_class_5genus18 = prune_taxa((tax_table(physeq_class18)[, "Genus.x"] %in% top5genus_nmds18), physeq_class18)



#REMOVING NA(s) FROM EACH LEVEL ####

## No Kingdom ####

## Phylum ####
physeq_class_5phyla18 = prune_taxa((tax_table(physeq_class18)[, "Phylum"] %in% top5phyla_nmds18), physeq_class18)
physeq_phylum5_no_na18 <- subset_samples(physeq_class_5phyla18, "Phylum" != "NA", )
data17.ord_phylum5_no_na18 <- ordinate(physeq_phylum5_no_na18, "NMDS", "bray")


## Class ####
physeq_class_5class18 = prune_taxa((tax_table(physeq_class18)[, "Class"] %in% top5class_nmds18), physeq_class18)
physeq_class5_no_na18 <- subset_samples(physeq_class_5class18, "Class" != "NA", )
data17.ord_class5_no_na18 <- ordinate(physeq_class5_no_na18, "NMDS", "bray")


## Order ####
physeq_class_5order18 = prune_taxa((tax_table(physeq_class18)[, "Order"] %in% top5order_nmds18), physeq_class18)
physeq_order5_no_na18 <- subset_samples(physeq_class_5order18, "Order" != "NA", )
data17.ord_order5_no_na18 <- ordinate(physeq_order5_no_na18, "NMDS", "bray")


## Family ####
physeq_class_5family18 = prune_taxa((tax_table(physeq_class18)[, "Family"] %in% top5family_nmds18), physeq_class18)
physeq_family5_no_na18 <- subset_samples(physeq_class_5family18, "Family" != "NA", )
data17.ord_family5_no_na18 <- ordinate(physeq_family5_no_na18, "NMDS", "bray")


## Genus ####
physeq_class_5genus18 = prune_taxa((tax_table(physeq_class18)[, "Genus.x"] %in% top5genus_nmds18), physeq_class18)
physeq_genus5_no_na18 <- subset_samples(physeq_class_5genus18, "Genus.x" != "NA", )
data17.ord_genus5_no_na18 <- ordinate(physeq_genus5_no_na18, "NMDS", "bray")





#NMDS of Top 5 Taxa w/out NAs ####

##Kingdom 
plot_ordination(physeq_class18, data18.ord, type="taxa", color="Kingdom")+
  scale_colour_manual(values=c("#FFA106","#2E86C1","#8FC172","#FF06F7", "#8D0101"))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS Kingdom Plots",
       subtitle = "Plotting Kingdom Diversity (Top 5)",
       caption = "Data source: Oyster 16s 2018")
ggsave(filename = "NMDS_Orig_King18.jpeg", plot=last_plot(), path ="Plot Diagrams 2018/Original_2018/NMDS /", width = 7, height = 5)  


#Phylum
plot_ordination(physeq_phylum5_no_na18, data17.ord_phylum5_no_na18, type="taxa", color="Phylum")+
  scale_colour_manual(values=c("#FFA106","#2E86C1","#8FC172","#FF06F7", "#8D0101"))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS Phylum Plots",
       subtitle = "Plotting Phylum Diversity (Top 5)",
       caption = "Data source: Oyster 16s 2018")
ggsave(filename = "NMDS_Orig_Phylum18.jpeg", plot=last_plot(), path ="Plot Diagrams 2018/Original_2018/NMDS /", width = 7, height = 5)  


#Class
plot_ordination(physeq_class5_no_na18, data17.ord_class5_no_na18, type="taxa", color="Class")+
  scale_colour_manual(values=c("#FFA106","#2E86C1","#8FC172","#FF06F7", "#8D0101"))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS Class Plots",
       subtitle = "Plotting Class Diversity (Top 5)",
       caption = "Data source: Oyster 16s 2018")
ggsave(filename = "NMDS_Orig_Class18.jpeg", plot=last_plot(), path ="Plot Diagrams 2018/Original_2018/NMDS /", width = 7, height = 5)  


#Order
plot_ordination(physeq_order5_no_na18, data17.ord_order5_no_na18, type="taxa", color="Order")+
  scale_colour_manual(values=c("#FFA106","#2E86C1","#8FC172","#FF06F7", "#8D0101"))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS Order Plots",
       subtitle = "Plotting Order Diversity (Top 5)",
       caption = "Data source: Oyster 16s 2018")
ggsave(filename = "NMDS_Orig_Order18.jpeg", plot=last_plot(), path ="Plot Diagrams 2018/Original_2018/NMDS /", width = 7, height = 5)  


#Family
plot_ordination(physeq_family5_no_na18, data17.ord_family5_no_na18, type="taxa", color="Family")+
  scale_colour_manual(values=c("#FFA106","#2E86C1","#8FC172","#FF06F7", "#8D0101"))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS Family Plots",
       subtitle = "Plotting Family Diversity (Top 5)",
       caption = "Data source: Oyster 16s 2018")
ggsave(filename = "NMDS_Orig_Family18.jpeg", plot=last_plot(), path ="Plot Diagrams 2018/Original_2018/NMDS /", width = 7, height = 5)  


#Genus
plot_ordination(physeq_genus5_no_na18, data17.ord_genus5_no_na18, type="taxa", color="Genus.x")+
  scale_colour_manual(values=c("#FFA106","#2E86C1","#8FC172","#FF06F7", "#8D0101"))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS Genus Plots",
       subtitle = "Plotting Genus Diversity (Top 5)",
       caption = "Data source: Oyster 16s 2018")
ggsave(filename = "NMDS_Orig_Genus18.jpeg", plot=last_plot(), path ="Plot Diagrams 2018/Original_2018/NMDS /", width = 7, height = 5)  




#NMDS of RTFM Significant ####
##physeq_sig18 (Phyloseq Object)####
##rftmdata18.ord (Ordination Object)####

#Isolating for the top 5 of significants
##Kingdom ####
kingdom.sig18 = tapply(taxa_sums(physeq_sig18), tax_table(physeq_sig18)[, "Kingdom"], sum, na.rm=TRUE)
top5kingdom_nmdsSIG18 = names(sort(kingdom.sig18, TRUE))[1:5]
physeq_SIG5kingdom18 = prune_taxa((tax_table(physeq_sig18)[, "Kingdom"] %in% top5kingdom_nmdsSIG18), physeq_sig18)

physeq_SigKing5_no_na18 <- subset_samples(physeq_SIG5kingdom18, "Kingdom" != "NA", )
data17.ord_SigKing5_no_na18 <- ordinate(physeq_SigKing5_no_na18, "NMDS", "bray")


##Phylum #### 
phylum.sig18 = tapply(taxa_sums(physeq_sig18), tax_table(physeq_sig18)[, "Phylum"], sum, na.rm=TRUE)
top5phyla_nmdsSIG18 = names(sort(phylum.sig18, TRUE))[1:5]
physeq_SIG5phyla18 = prune_taxa((tax_table(physeq_sig18)[, "Phylum"] %in% top5phyla_nmdsSIG18), physeq_sig18)

physeq_SigPhylum5_no_na18 <- subset_samples(physeq_SIG5phyla18, "Phylum" != "NA", )
data17.ord_SigPhylum5_no_na18 <- ordinate(physeq_SigPhylum5_no_na18, "NMDS", "bray")


##Class ####
class.sig18 = tapply(taxa_sums(physeq_sig18), tax_table(physeq_sig18)[, "Class"], sum, na.rm=TRUE)
top5class_nmdsSIG18 = names(sort(class.sig18, TRUE))[1:5]
physeq_SIG5class18 = prune_taxa((tax_table(physeq_sig18)[, "Class"] %in% top5class_nmdsSIG18), physeq_sig18)

physeq_SigClass5_no_na18 <- subset_samples(physeq_SIG5class18, "Class" != "NA", )
data17.ord_SigClass5_no_na18 <- ordinate(physeq_SigClass5_no_na18, "NMDS", "bray")


##Order #### 
order.sig18 = tapply(taxa_sums(physeq_sig18), tax_table(physeq_sig18)[, "Order"], sum, na.rm=TRUE)
top5order_nmdsSIG18 = names(sort(order.sig18, TRUE))[1:5]
physeq_SIG5order18 = prune_taxa((tax_table(physeq_sig18)[, "Order"] %in% top5order_nmdsSIG18), physeq_sig18)

physeq_SigOrder5_no_na18 <- subset_samples(physeq_SIG5order18, "Order" != "NA", )
data17.ord_SigOrder5_no_na18 <- ordinate(physeq_SigOrder5_no_na18, "NMDS", "bray")


##Family ####
family.sig18 = tapply(taxa_sums(physeq_sig18), tax_table(physeq_sig18)[, "Family"], sum, na.rm=TRUE)
top5family_nmdsSIG18 = names(sort(family.sig18, TRUE))[1:5]
physeq_SIG5family18 = prune_taxa((tax_table(physeq_sig18)[, "Family"] %in% top5family_nmdsSIG18), physeq_sig18)

physeq_SigFamily5_no_na18 <- subset_samples(physeq_SIG5family18, "Order" != "NA", )
data17.ord_SigFamily5_no_na18 <- ordinate(physeq_SigFamily5_no_na18, "NMDS", "bray")


##Genus ####
genus.sig18 = tapply(taxa_sums(physeq_sig18), tax_table(physeq_sig18)[, "Genus.x"], sum, na.rm=TRUE)
top5genus_nmdsSIG18 = names(sort(genus.sig18, TRUE))[1:5]
physeq_SIG5genus18 = prune_taxa((tax_table(physeq_sig18)[, "Genus.x"] %in% top5genus_nmdsSIG18), physeq_sig18)

physeq_SigGenus5_no_na18 <- subset_samples(physeq_SIG5genus18, "Genus.x" != "NA", )
data17.ord_SigGenus5_no_na18 <- ordinate(physeq_SigGenus5_no_na18, "NMDS", "bray")


#Kingdom 
physeq_SigKing5_no_na18 <- subset_samples(physeq_SIG5kingdom18, "Kingdom" != "NA", )
data17.ord_SigKing5_no_na18 <- ordinate(physeq_SigKing5_no_na18, "NMDS", "bray")

plot_ordination(physeq_SigKing5_no_na18, data17.ord_SigKing5_no_na18, type="taxa", color="Kingdom")+
  scale_colour_manual(values=c("#FFA106","#2E86C1","#8FC172","#FF06F7", "#8D0101"))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#2E86C1"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS Kingdom Plots",
       subtitle = "Plotting Kingdom Diversity (RTFM Significant OTUs)",
       caption = "Data source: Oyster 16s 2018")
ggsave(filename = "NMDS_RFTMSig_King18.jpeg", plot=last_plot(), path ="Plot Diagrams 2018/RFTM Significant/", width = 7, height = 5)  


#Phylum
physeq_SigPhylum5_no_na18 <- subset_samples(physeq_SIG5phyla18, "Phylum" != "NA", )
data17.ord_SigPhylum5_no_na18 <- ordinate(physeq_SigPhylum5_no_na18, "NMDS", "bray")

plot_ordination(physeq_SigPhylum5_no_na18, data17.ord_SigPhylum5_no_na18, type="taxa", color="Phylum")+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#2E86C1"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS Phylum Plots",
       subtitle = "Plotting Phylum Diversity (RTFM Significant OTUs)",
       caption = "Data source: Oyster 16s 2018")
ggsave(filename = "NMDS_RFTMSig_Phylum18.jpeg", plot=last_plot(), path ="Plot Diagrams 2018/RFTM Significant/", width = 7, height = 5)  


#Class
physeq_SigClass5_no_na18 <- subset_samples(physeq_SIG5class18, "Class" != "NA", )
data17.ord_SigClass5_no_na18 <- ordinate(physeq_SigClass5_no_na18, "NMDS", "bray")

plot_ordination(physeq_SigClass5_no_na18, data17.ord_SigClass5_no_na18, type="taxa", color="Class")+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#2E86C1"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS Class Plots",
       subtitle = "Plotting Class Diversity (RTFM Significant OTUs)",
       caption = "Data source: Oyster 16s 2018")
ggsave(filename = "NMDS_RFTMSig_Class18.jpeg", plot=last_plot(), path ="Plot Diagrams 2018/RFTM Significant/", width = 7, height = 5)  


#Order
physeq_SigOrder5_no_na18 <- subset_samples(physeq_SIG5order18, "Order" != "NA", )
data17.ord_SigOrder5_no_na18 <- ordinate(physeq_SigOrder5_no_na18, "NMDS", "bray")

plot_ordination(physeq_SigOrder5_no_na18, data17.ord_SigOrder5_no_na18, type="taxa", color="Order")+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#2E86C1"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS Order Plots",
       subtitle = "Plotting Order Diversity (RTFM Significant OTUs)",
       caption = "Data source: Oyster 16s 2018")
ggsave(filename = "NMDS_RFTMSig_Order18.jpeg", plot=last_plot(), path ="Plot Diagrams 2018/RFTM Significant/", width = 7, height = 5)  


#Family
physeq_SigFamily5_no_na18 <- subset_samples(physeq_SIG5family18, "Order" != "NA", )
data17.ord_SigFamily5_no_na18 <- ordinate(physeq_SigFamily5_no_na18, "NMDS", "bray")

plot_ordination(physeq_SigFamily5_no_na18, data17.ord_SigFamily5_no_na18, type="taxa", color="Family")+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#2E86C1"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS Family Plots",
       subtitle = "Plotting Family Diversity (RTFM Significant OTUs)",
       caption = "Data source: Oyster 16s 2018")
ggsave(filename = "NMDS_RFTMSig_Family18.jpeg", plot=last_plot(), path ="Plot Diagrams 2018/RFTM Significant/", width = 7, height = 5)  


#Genus - top 5 significant 
physeq_SigGenus5_no_na18 <- subset_samples(physeq_SIG5genus18, "Genus.x" != "NA", )
data17.ord_SigGenus5_no_na18 <- ordinate(physeq_SigGenus5_no_na18, "NMDS", "bray")

plot_ordination(physeq_SigGenus5_no_na18, data17.ord_SigGenus5_no_na18, type="taxa", color="Genus.x")+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#2E86C1"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS Genus Plots",
       subtitle = "Plotting Genus Diversity (RTFM Significant OTUs)",
       caption = "Data source: Oyster 16s 2018")
ggsave(filename = "NMDS_RFTMSig_Genus18.jpeg", plot=last_plot(), path ="Plot Diagrams 2018/RFTM Significant/", width = 7, height = 5)  


#Genus - Significant 
plot_ordination(physeq_sig18, rftmdata18.ord, type="taxa", color="Genus.x")+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#2E86C1"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "NMDS Genus Plots",
       subtitle = "Plotting Genus Diversity (RTFM Significant OTUs)",
       caption = "Data source: Oyster 16s 2018")
ggsave(filename = "NMDS_RFTMSig_Genus18.jpeg", plot=last_plot(), path ="Plot Diagrams 2018/RFTM Significant/", width = 7, height = 5)  









#Evaluation of Samples ####

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




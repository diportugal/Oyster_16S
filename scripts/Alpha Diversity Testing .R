#Graphical Analysis - Testing Alpha Diversity 
#2021-06-24
#Author:Diana Portugal
#Contact: dportugal8@gmail.com 


## LOADING REQUIRED R PACKAGES ####
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("plyr"); packageVersion("plyr")
library("data.table")
theme_set(theme_bw())
data("GlobalPatterns")


####
#Notes:
#Alpha diversity is used to measure the diversity within a sample. "How many?"
#High Alpha diversity = many organisms 
#Considers species richness and species evennes



physeq_class = phyloseq(OTU, TAX, SAMP) 
physeq_class


physeq_class <- prune_taxa(taxa_sums(GlobalPatterns) > 0, GlobalPatterns)

plot_richness(physeq_class)
#This brings up the 7 plot types: Observed, Chao1, ACE, Shannon, Simpson, InvSimpson, and Fisher

plot_richness(physeq_class, measures=c("Chao1", "Shannon"))
#This brings up only the plot type that you want out of the six

plot_richness(physeq_class, x="samples", measures=c("Chao1", "Shannon"))
#Taking the samples on the x axis



sample_data(physeq_class)$Site.x <- get_variable(physeq_class, "Site.x") 

plot_richness(physeq_class, x="Site.x", color="Site.x", measures=c("Chao1", "Shannon"))



Phy.ord <- ordinate(physeq_class, "NMDS", "bray")

plot_richness(physeq_class, Phy.ord, x="samples", color="samples", measures=c("Chao1", "Shannon"))



??plot_richness()





#Use "require(RColorBrewer)"


#Sarah's notes:
#plot_bar(metaSF,  fill="Category", x="Replicate") +
  geom_bar(aes(color=Category, fill=Category), stat="identity", position="stack")+
  facet_grid(Year~Site_Name, scales="free_x")+
  scale_fill_manual(values=mycolors)+
  scale_color_manual(values=mycolors)+
  theme_bw()+
  theme(legend.position="top", legend.text=element_text(size=10), panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor.x = element_blank(),
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10))







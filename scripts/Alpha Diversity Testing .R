p = plot_bar(ent10, "Genus", fill="Genus", facet_grid=SeqTech~Enterotype)
p + geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")#Graphical Analysis - Testing Alpha Diversity 
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
#Considers species richness and species evenness



physeq_class = phyloseq(OTU, TAX, SAMP) 
physeq_class

pruned_phydata17 <- prune_taxa(taxa_sums(GlobalPatterns) > 0, GlobalPatterns)


plot_richness(pruned_phydata17)
#This brings up the 7 plot types: Observed, Chao1, ACE, Shannon, Simpson, InvSimpson, and Fisher

plot_richness(pruned_phydata17, measures=c("Chao1", "Shannon"))
#This brings up only the plot type that you want out of the six

plot_richness(pruned_phydata17, x="samples", measures=c("Chao1", "Shannon"), title = "Testing Testing")
#Taking the samples on the x axis

plot_richness(OTU,taxa_are_rows(physeq), x="samples", color="Treatment2", measures=c("Chao1", "Shannon"))




require(RColorBrewer)


??plot_richness()
??taxa_are_rows


plot_bar(pruned_phydata17,  fill="Treatment2", x="Replicate") +
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







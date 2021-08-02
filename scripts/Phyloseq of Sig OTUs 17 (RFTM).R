#Oyster 16S - Phyloseq Anaysis for 2017 RFTM Significant Data
#Author: Diana Portugal 
#Contact: dportugal8@gmail.com 
#Date 24-7-2021



library("phyloseq"); packageVersion("phyloseq")
data("GlobalPatterns")
library("ggplot2"); packageVersion("ggplot2")
library("plyr"); packageVersion("plyr")
theme_set(theme_bw())  shagdhagsd
library("data.table")
library("vegan")



## Isolating Significant OTUs for RFTM Variable (2017) ####

#sig_OTUs_RFTM17 <- RFTM_sig17 %>% 
  #select(Kingdom, Phylum, Class, Order, Family, Genus.x)

rftmpos = RFTM_sig17[RFTM_sig17$log2FoldChange<0,]
rftmasv = row.names(rftmpos)

phylo17_rftm = prune_taxa(taxa_sums(physeq_class) > 0,physeq_class ) 

rftmdata17.ord <- ordinate(phylo17_rftm, "NMDS", "bray")



#ALPHA DIVERSITY ####
plot_richness(phylo17_rftm)
#What does plot richness take in as the variable it is evaluating?

plot_richness(phylo17_rftm, measures=c("Chao1", "Shannon"))
#Selecting these two methods of evaluating 

p <- plot_richness(phylo17_rftm, x="RFTM_score.x", measures=c("Chao1", "Shannon"), col="Site.x")+
  geom_point(size=2, alpha=0.7)
p
#selecting the variable "RFTM_score.x" and adding a larger point with transparency 

p$layers <- p$layers[-1]
#removing 

p+
  geom_point(size = 5, alpha = 0.7)+
  stat_summary(fun.data = mean_se, geom ="errorbar", width=0.2)



dim(phylo17_rftm)


#I wanted to compare the original Phyloseq object to the one with the significant OTUs from RFTM DESeq
#They appear to be the same?? 
plot_richness(phylo17_rftm, x="RFTM_score.x", measures=c("Simpson", "Shannon"), col="Site.x")+
  geom_point(size=2, alpha=0.7)+ 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5,))+
  labs(title = "Diversity of RFTM Scores OTUs - Significant OTUs",
       subtitle = "Measuring alpha diversity of RFTM score with Shannon and Simpson",
       caption = "Data source: Oyster 16s 2017")
ggsave("Diversity_RFTM_Sig.jpeg", path="Plot Diagrams/", width = 7, height = 5)


plot_richness(physeq_class, x="RFTM_score.x", measures=c("Simpson", "Shannon"), col="Site.x")+
  geom_point(size=2, alpha=0.7)+ 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5,))+
  labs(title = "Diversity of RFTM Scores - Original Data",
       subtitle = "Measuring alpha diversity of RFTM score with Shannon and Simpson",
       caption = "Data source: Oyster 16s 2017")
ggsave("Diversity_RFTM_Orig.jpeg", path="Plot Diagrams/", width = 7, height = 5)


plot_richness(phylo17_rftm, x="RFTM_score.x", measures=c("Simpson", "Shannon"), col="Site.x")+
  geom_point(size=2, alpha=0.7)+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5,))+
  labs(title = "Diversity of RFTM Scores by OTUs",
       subtitle = "Measuring alpha diversity of RFTM score with Shannon and Simpson",
       caption = "Data source: Oyster 16s 2017")






plot_bar(physeq_class, aes(x="Site.x", fill="Species.x"))

plot_bar(physeq_class,"Site.x")


plot_bar(physeq_class,"Site.x", fill="Species.x")


plot_bar(physeq_class,"RFTM_score.x", fill="Site.x")+
  geom_col()


plot_bar(physeq_class,"RFTM_score.x", fill="Site.x")+
  geom_col()+
  geom_text(aes(label="Site.x"), position=position_stack(vjust=0.5), colour="white", size=5)

plot_bar(physeq_class,"RFTM_score.x", fill="Site.x")+
  geom_col()




parse_phyloseq








































































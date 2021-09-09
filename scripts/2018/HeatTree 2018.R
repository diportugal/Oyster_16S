#HEAT TREE
#2021-08-02
#Author:Diana Portugal
#Contact: dportugal8@gmail.com 

#install.packages("metacoder")
library(metacoder)
library(phyloseq)


#Heat Tree with Original Data ####

##Creating a taxmap object from original Phyloseq object ####
taxmap18 = parse_phyloseq(physeq_class18)

##Creating a taxmap object from significant Phyloseq object ####
taxmapsig18 <- parse_phyloseq(physeq_sig18)



#Generating a heat tree
taxmapsig18 %>% 
  heat_tree(node_label = taxon_names,
            node_size = n_obs(taxmapsig18),
            node_color = n_obs(taxmapsig18),
            layout = "da", initial_layout = "re",
            title = "RFTM Microbiome Taxonomy 2018 Oyster 16S")
ggsave(filename = "HeatTree_Significant_18.jpeg", plot=last_plot(), path="Plot Diagrams 2018/RFTM Significant/", width = 7, height = 5)  
#figure out how to make edits to the text and color 


taxmapsig18 %>% 
  heat_tree(node_label = taxon_names,
            node_size = n_obs(taxmapsig18),
            node_color = n_obs(taxmapsig18),
            node_color_range = diverging_palette(),
            layout = "da", initial_layout = "re",
            title = "RFTM Microbiome Taxonomy 2018 Oyster 16S XXX")


#comment 

#HEAT TREE
#2021-07-30
#Author:Diana Portugal
#Contact: dportugal8@gmail.com 

#install.packages("metacoder")
library(metacoder)
library(phyloseq)



#Creating a taxmap object from original Phyloseq object ####
taxmap <- parse_phyloseq(physeq_class)

#Creating a taxmap object from significant Phyloseq object (RFTM_score.x) ####
taxmapsig <- parse_phyloseq(phylo17_rftm)

#Creating a taxmap object from significant Phyloseq object (PeaCrab) ####
taxmapPeaSig <- parse_phyloseq(physeq_sigPea)



## Simple Heat Trees (PeaCrab) ####
taxmapPeaSig %>% 
  heat_tree(node_label = taxon_names,
            node_size = n_obs(taxmapPeaSig),
            node_color = n_obs(taxmapPeaSig),
            layout = "da", initial_layout = "re",
            title = "Microbiome Taxonomy 2017 Oyster 16S")
ggsave(filename = "HeatTree_PeaSig_17.jpeg", plot=last_plot(), path="Plot Diagrams/Peacrab/", width = 7, height = 5)  


## Simple Heat Trees (Original) ####
taxmap %>% 
  heat_tree(node_label = taxon_names,
            node_size = n_obs(taxmap),
            node_color = n_obs(taxmap),
            layout = "da", initial_layout = "re")
ggsave(filename = "HeatTree_Original_17.jpeg", plot=last_plot(), path="Plot Diagrams/NMDS Original /", width = 7, height = 5)  
#figure out how to make edits to the text and color 


## Simple Heat Trees (Significant) ####
taxmapsig %>% 
  heat_tree(node_label = taxon_names,
            node_size = n_obs(taxmapsig),
            node_color = n_obs(taxmapsig),
            layout = "da", initial_layout = "re")
ggsave(filename = "HeatTree_RFTMSIG_17.jpeg", plot=last_plot(), path="Plot Diagrams/", width = 7, height = 5)  




#Trying out a different layout style (Dont Use) #### 
taxmap %>% 
  heat_tree(node_label = taxon_names,
            node_size = n_obs(taxmap),
            node_color = n_obs(taxmap),
            layout = "davidson-harel", initial_layout = "reingold-tilford",
            title = "Microbiome Taxonomy 2017 Oyster 16S")
ggsave(filename = "HeatTree_Original2_17.jpeg", plot=last_plot(), path="Plot Diagrams/", width = 7, height = 5)  
#This one is a bit more crowded.



#Trying out a different Color style (Dont Use) #### 
taxmapsig %>% 
  heat_tree(node_label = taxon_names,
            node_size = n_obs(taxmapsig),
            node_color = n_obs(taxmapsig),
            node_color_range = diverging_palette(),
            layout = "da", initial_layout = "re",
            title = "Microbiome Taxonomy 2017 Oyster 16S (Sig)")
ggsave(filename = "HeatTree_RFTMSIG_17.jpeg", plot=last_plot(), path="Plot Diagrams/", width = 7, height = 5)  





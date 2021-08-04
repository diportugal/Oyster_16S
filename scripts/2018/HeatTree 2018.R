#HEAT TREE
#2021-08-02
#Author:Diana Portugal
#Contact: dportugal8@gmail.com 

#install.packages("metacoder")
library(metacoder)
library(phyloseq)


#Heat Tree with Original Data ####

##Creating a taxmap object from original Phyloseq object ####
taxmap18 <- parse_phyloseq(physeq_class18)

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



#Trying out a different layout style 
taxmap %>% 
  heat_tree(node_label = taxon_names,
            node_size = n_obs(taxmap),
            node_color = n_obs(taxmap),
            layout = "davidson-harel", initial_layout = "reingold-tilford",
            title = "Microbiome Taxonomy 2017 Oyster 16S")
ggsave(filename = "HeatTree_Original2_17.jpeg", plot=last_plot(), path="Plot Diagrams/", width = 7, height = 5)  
#This one is a bit more crowded.



#Heat Tree with Significant Data ####
taxmapsig %>% 
  heat_tree(node_label = taxon_names,
            node_size = n_obs(taxmapsig),
            node_color = n_obs(taxmapsig),
            layout = "da", initial_layout = "re",
            title = "Microbiome Taxonomy 2017 Oyster 16S (Sig)")
ggsave(filename = "HeatTree_RFTMSIG_17.jpeg", plot=last_plot(), path="Plot Diagrams/", width = 7, height = 5)  


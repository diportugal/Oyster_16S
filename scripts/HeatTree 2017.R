#HEAT TREE
#2021-07-30
#Author:Diana Portugal
#Contact: dportugal8@gmail.com 

#install.packages("metacoder")
library(metacoder)
library(phyloseq)



#Creating a taxmap object from original Phyloseq object ####
taxmap <- parse_phyloseq(physeq_class)
physeq_class
tax_table(physeq_class)

#Creating a taxmap object from significant Phyloseq object (RFTM_score.x) ####
taxmapRFTMsig <- parse_phyloseq(phylo17_rftm)

phylo17_rftm
tax_table(phylo17_rftm)

taxmapRFTMsig$data$deseq=RFTM_sig17

#Creating a taxmap object from significant Phyloseq object (PeaCrab) ####
taxmapPeaSig <- parse_phyloseq(phylo17_pea)
tax_table(phylo17_pea)
phylo17_pea


## Simple Heat Trees (PeaCrab) ####
taxmapPeaSig %>% 
  heat_tree(node_label = taxon_names,
            node_size = n_obs(taxmapPeaSig),
            node_color = n_obs(taxmapPeaSig),
            layout = "da", initial_layout = "re")+
  ggtitle("OYSTER - PEACRAB HEAT TREE")+
  theme(plot.title=element_text(hjust=0.5, vjust=0.7, face="bold", colour='Black', size=10))
ggsave(filename = "HT_PeaSig_17.jpeg", plot=last_plot(), path="Plot Diagrams/", width = 7, height = 5)  



## Simple Heat Trees (Original) ####
taxmap %>% 
  heat_tree(node_label = taxon_names,
            node_size = n_obs(taxmap),
            node_color = n_obs(taxmap),
            layout = "da", initial_layout = "re")+
  ggtitle("OVERALL DATA HEAT TREE")+
  theme(plot.title=element_text(hjust=0.5, vjust=0.7, face="bold", colour='Black', size=10))
ggsave(filename = "HT_Original_17.jpeg", plot=last_plot(), path="Plot Diagrams/", width = 7, height = 5)  


## Simple Heat Trees (RFTM Sig) ####
taxmapRFTMsig %>% 
  heat_tree(node_label = taxon_names,
            node_size = n_obs(taxmapRFTMsig),
            node_color = n_obs(taxmapRFTMsig),
            layout = "da", initial_layout = "re")+
  ggtitle("OYSTER - P. MARINUS HEAT TREE")+
  theme(plot.title=element_text(hjust=0.5, vjust=0.5, face="bold", colour='Black', size=10))
ggsave(filename = "HeatTree_RFTMSIG_17.jpeg", plot=last_plot(), path="Plot Diagrams/", width = 7, height = 5)  




#Trying out a different layout style (Dont Use) #### 
taxmapRFTMsig %>% 
  heat_tree(node_label = taxon_names,
            node_size = n_obs(taxmapRFTMsig),
            node_color = n_obs(taxmapRFTMsig),
            layout = "davidson-harel", initial_layout = "reingold-tilford")
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





#HEAT TREE

#install.packages("metacoder")
library(metacoder)
library(phyloseq)

#Heat Tree

physeq_class %>% 
  heat_tree(node_label = gsub(pattern = "\\[|\\]", replacement = "", taxon_names),
            node_size = ,
            node_color = ,
            node_color_axis_label = "OTU count",
            layout = "davidson-harel", initial_layout = "reingold-tilford")



heat_tree(physeq_class)
heat_tree(TAX$Kingdom)


taxmap <- parse_phyloseq(physeq_class)

taxmap %>% 
  heat_tree(node_label = "taxon_names",
            node_size = Class,
            node_color = Class,
            layout = "da", initial_layout = "re",
            title = "Taxa in Oysters")



(node_label = taxon_names,
  node_size = leaf,
  node_color = leaf,
  layout = "da", initial_layout = "re",
  title = "Taxa in leafs")








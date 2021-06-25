#Graphical Analysis - Testing 
#2021-06-23
#Author:Diana Portugal
#Contact: dportugal8@gmail.com 


## LOADING REQUIRED R PACKAGES ####
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("plyr"); packageVersion("plyr")
library("data.table")
theme_set(theme_bw())



## MAKING A PLOT GRAPH ####
plot_bar(physeq_class, x= "Treatment2", y="RFTM_score.x", fill = "Species.x")
#NOTE - The graphs take a while to load 



## MAKING AN NMDS PLOT (SCATTED DOTS): JUST OTUs ####
Phy.ord <- ordinate(physeq_class, "NMDS", "bray")
plot1 = plot_ordination(physeq_class, Phy.ord, type="taxa", color="Kingdom", title="TAXA")
print(plot1)

Phy.ord <- ordinate(physeq_class, "NMDS", "bray")
plot1 = plot_ordination(physeq_class, Phy.ord, type="taxa", color="Kingdom", title="TAXA")
print(plot1)

#Grouping by Phylum and doing 3 plots stacked on top of each other
plot1 + facet_wrap(~Kingdom, 3) 

#QUESTION: CAN WE CHANGE THE COLORS?!?!?



## MAKING A PLOT: JUST SAMPLES ####
plot2 = plot_ordination(physeq_class, Phy.ord, type="samples", color= sample_variables(physeq), shape=NULL) 
plot2 + geom_polygon(aes(fill=NULL)) + geom_point(size=2, color= "red") + ggtitle("PLOT 2: JUST SAMPLES")
#LOOK MORE INTO MANIPULATING THESE VARIABLES



## MAKING A BIPLOT GRAPHIC ####
plot3 = plot_ordination(physeq_class, Phy.ord, type="biplot", color=NULL, shape=NULL, title="PLOT 3: BIPLOT GRAPHIC")
# Some stuff to modify the automatic shape scale, not nessary to include but look into!!
GP1.shape.names = get_taxa_unique(GP1, "Phylum")
GP1.shape <- 15:(15 + length(GP1.shape.names) - 1)
names(GP1.shape) <- GP1.shape.names
GP1.shape["samples"] <- 16
p3 + scale_shape_manual(values=GP1.shape)

print(plot3)


## MAKING A SPLIT GRAPHIC ####

plot4 = plot_ordination(physeq_class, Phy.ord, type="split", color=NULL, shape=NULL, label=NULL, title="PLOT 4 : SPLIT GRAPHIC") 
print(plot4)


??plot_ordination()











#ANSWERING QUESTIONS

#plot_q1 <- plot_bar(physeq_class, x= "Treatment2", y="RFTM_score.x", fill = "Species.x")

phylum.sum = tapply(taxa_sums(physeq_class), tax_table(physeq_class)[, "Phylum"], sum, na.rm=TRUE)
top5phyla = names(sort(phylum.sum, TRUE))[1:5]
physeq_class = prune_taxa((tax_table(physeq_class)[, "Phylum"] %in% top5phyla), physeq_class)


Phy.ord <- ordinate(physeq_class, "NMDS", "bray")
plot_q1 = plot_ordination(physeq_class, Phy.ord, type="taxa", color="Phylum", title="Plot Q1")
print(plot_q1)


plot2_q1 = plot_ordination(physeq_class, Phy.ord, type="samples", color="Treatment2", shape="Site.y") 
plot2_q1 + geom_polygon(aes(fill=Treatment2)) + geom_point(size=5) + ggtitle("Treatment2 x Site")


#####

#Site x treatment 2 does not show any differences
#Site x and Site y are the same

#####


plot2_q1 = plot_ordination(physeq_class, Phy.ord, type="samples", color="Treatment2", shape="Site.y") 
plot2_q1 + geom_polygon(aes(fill=Treatment2)) + geom_point(size=10) + ggtitle("samples")


plot3_q1 = plot_ordination(physeq_class, Phy.ord, type="samples", color="Treatment2", shape="Site.y") 
plot3_q1 + geom_polygon(aes(fill=Treatment2)) + geom_point(size=10) + ggtitle("samples")

plot3_q1 + facet_wrap(~Phylum, 3)


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
data("enterotype")



# MAKING A PLOT GRAPH ####
plot_bar(physeq_class, x= "Treatment2", y="Weight_post", fill = "#FF9BD4")
#NOTE - The graphs take a while to load 
#NOTE - The y axis is in a different order of magnitude




# PLOT ORDINATION ####

#THESE ARE NMDS CHARTS
#NMDS CHART DISPLAYING DIFFERENCES BY TREATMENT2 TYPES
Phy.ord <- ordinate(physeq_class, "NMDS", "bray")
nmds_t2 = plot_ordination(physeq_class, Phy.ord, type ="samples", 
                        color = "Treatment2", label = NULL, title = "NMDS GRAPH OF TREATMENT2") 
print(nmds_t2)


#NMDS CHART DISPLAYING TAXA DIFFERENCES BY KINGDOM
nmds_taxa = plot_ordination(physeq_class, Phy.ord, type="taxa", color="Kingdom", title="NMDS GRAPH OF TAXA")
print(nmds_taxa)




#FACET WRAPS 
#Grouping by Phylum and doing 3 plots stacked on top of each other
plot1 + facet_wrap(~Kingdom, 3) 



plot_samp = plot_ordination(physeq_class, Phy.ord, type="samples", color= "Treatment2", shape="Site.x") 
plot_samp


## MAKING A PLOT: JUST SAMPLES ####
#This plot is graphing Site/Treatment2 NMDS
plot2 = plot_ordination(physeq_class, Phy.ord, type="samples", color= "Treatment2", shape="Site.x") 
plot2 + geom_polygon(aes(fill=NULL)) + geom_point(size=2, color= "#33B5FF") + ggtitle("PLOT 2: Site/Treatment2 NMDS")





??plot_ordination()


## MAKING A BIPLOT GRAPHIC ####
plot3 = plot_ordination(physeq_class, Phy.ord, type="biplot", color="NULL", shape=NULL, title="PLOT 3: BIPLOT GRAPHIC")
print(plot3)

# Some stuff to modify the automatic shape scale, not nessary to include but look into!!
physeq_class.shape.names = get_taxa_unique(physeq_class, "Kingdom")
physeq_class.shape <- 15:(15 + length(physeq_class.shape.names) - 1)
names(physeq_class.shape) <- physeq_class.shape.names
physeq_class.shape["samples"] <- 16
plot3 + scale_shape_manual(values=physeq_class.shape)





#FILTERING FOR THE TOP 10 MOST COMON TAXA
TopNOTUs <- names(sort(taxa_sums(enterotype), TRUE)[1:10])
ent10   <- prune_taxa(TopNOTUs, enterotype)
#


p3 = plot_ordination(physeq_class, Phy.ord, type="biplot", color="Treatment2", shape="Species.x", title="P3: IM TRYING")
# Some stuff to modify the automatic shape scale, not nessary to include but look into!!
physeq_class.shape.names = get_taxa_unique(physeq_class, "Kingdom")
physeq_class.shape <- 15:(15 + length(physeq_class.shape.names) - 1)
names(physeq_class.shape) <- physeq_class.shape.names
physeq_class.shape["ent10"] <- 16
p3 + scale_shape_manual(values=physeq_class.shape)

print(p3)







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




## BAR PLOTS ####

data("physeq_class")
physeq_class = subset_taxa(physeq_class ) 

plot_bar(physeq_class)


plot_bar(physeq_class, fill="Phylum")

plot_bar(physeq_class, x="Site.x", fill="Phylum")
#Why does it not want to take the site in which the sample was taken



#Histograms with basic R

hist(metadata17_df$RFTM_score.x, )

hist(metadata17_df$peacrabs.x, breaks ="sturges")


hist(metadata17_df$delta_weight17) 



??hist()












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





#SIGNIFICANT OTUs from DESeq2 ####

hist()


#Histogram of Significant OTUs from the DESeq Analysis 
hist(SigOTUs_RFTM_Tax17, 
     main="Significant OTUs from RFTM/Taxa",
     xlab="RFTM Score",
     col="darkmagenta")+
  geom_histogram(aes(fill = ))




plotDispEsts(RFTM_dds17)


plotMA( RFTM_res17, ylim = c(-1, 1) )

#The plot represents each gene with a dot. 
#The x axis is the average expression over all samples
#The y axis the log2 fold change between treatment and control
#Genes with an adjusted p value below a threshold (here 0.1, the default) are shown in red.



hist(RFTM_res17$pvalue, breaks=20, col="grey" )

#So you have a peak at 0, just like you saw in (A)… but you also have a peak close to 1. What do you do?
#Don’t apply false discovery rate control to these p-values yet. (Why not? Because some kinds of FDR control are based on the assumption that your p-values near 1 are uniform. If you break this assumption, you’ll get way fewer significant hypotheses. Everyone loses).
#Instead, figure out why your p-values show this behavior, and solve it appropriately:
#Are you applying a one-tailed test (for example, testing whether each gene increased its expression in response to a drug)? If so, those p-values close to 1 are cases that are significant in the opposite direction (cases where genes decreased their expression). If you want your test to find these cases, switch to a two-sided test. If you don’t want to include them at all, you can try filtering out all cases where your estimate is in that direction.
#Do all the p-values close to 1 belong to some pathological case? An example from my own field: RNA-Seq data, which consists of read counts per gene in each a variety of conditions, will sometimes include genes for which there are no reads in any condition. Some differential expression software will report a p-value of 1 for these genes. If you can find problematic cases like these, just filter them out beforehand (it’s not like you’re losing any information!)





hist(RFTM_res17$padj, breaks=20, col="grey" )




RFTM_sig17 #Sig OTUs 




#make New Physeq object with significant  OTU 

physeq_RFTM_Sig = phyloseq(RFTM_sig17)


#Merge with Taxa 

#Prune taxa to keep only Taxa from OTUs

#The OTUs will prune, you will only keep the OTUs in the taxonomy table

#prune taxonomy by sig otus




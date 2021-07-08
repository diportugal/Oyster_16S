#Oyster 16S Phyloseq Barplot Anaysis for 2018 Data
#Author: Diana Portugal 
#Contact: dportugal8@gmail.com 

library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
theme_set(theme_bw())





plot_bar(physeq_class18, fill="Kingdom")

#This loooks cray 



data("physeq_class18")
pc18.sp = subset_taxa(physeq_class18, Phylum == "Firmicutes")

plot_bar(physeq_class_5genus18, x="Species.x", fill="Genus")


plot_bar(subset_taxa(physeq_class18, Phylum == "Verrucomicrobiota"))



#Kingdom per Species ID
plot_bar(physeq_class18, x="Species.x", fill="Kingdom")+
  geom_col()+
  facet_wrap(~Kingdom, )


plot_bar(physeq_class18, x="Species.x", fill="Kingdom")+       #This plot is far less clear than ^^^
  geom_col()+
  facet_wrap(~Species.x, )


#Top 5 Phyla per Species ID
plot_bar(physeq_class_5phyla18, x="Species.x", fill="Phylum")+
  geom_col()

  
#Top 5 Class per Species ID
plot_bar(physeq_class_5class18, x="Species.x", fill="Class")+
  geom_col()+
  scale_fill_manual(values=c("#008ae6", "#00b300", "#ff9900", "#ff4d4d", "#bf80ff"))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "Relative Abundance Microbes by Class",
       subtitle = "Top Five Class evaluated",
       caption = "Data source: Oyster 16s 2018")


  
#Top 5 Order per Species ID
plot_bar(physeq_class_5order18, x="Species.x", fill="Order")+
  geom_col()+
  scale_fill_manual(values=c("#008ae6", "#00b300", "#ff9900", "#ff4d4d", "#bf80ff"))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "Relative Abundance Microbes by Order",
       subtitle = "Top Five Order evaluated",
       caption = "Data source: Oyster 16s 2018")

  
#Top 5 Family per Species ID
plot_bar(physeq_class_5family18, x="Species.x", fill="Family")+
  geom_col()+
  scale_fill_manual(values=c("#008ae6", "#00b300", "#ff9900", "#ff4d4d", "#bf80ff"))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "Relative Abundance Microbes by Family",
       subtitle = "Top Five Family evaluated",
       caption = "Data source: Oyster 16s 2018")


  
#Top 5 Genus per Species ID
plot_bar(physeq_class_5genus18, x="Species.x", fill="Genus")+
  geom_col()+
  scale_fill_manual(values=c("#008ae6", "#00b300", "#ff9900", "#ff4d4d", "#bf80ff"))+
  theme(legend.position="right", legend.text=element_text(size=10), 
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10), 
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15, colour = "#4E84C4"), 
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "Relative Abundance Microbes by Genus",
       subtitle = "Top Five Genus evaluated",
       caption = "Data source: Oyster 16s 2018")



#Question:Fixing abundance scale?? ####

plot_bar(physeq_class18, x="Species.x", fill="RFTM_score.x")+
  geom_col()


plot_bar(physeq_class18, x="Species.x", fill="delta_weight18")+
  geom_col()





#Oyster 16S Metadata Cleaning 
#Author: Diana Portugal 
#Contact: dportugal8@gmail.com 

#KEY ####
#Renaming the species collected
#CV - Crassostrea virginica (Oyster)
#IR - Ischadium recurvum (‘hooked’ Mussel)
#MM - Macoma mitchelli (Clam) NEW NAME is AM - Ameritella mitchelli
#MB = Macoma balthica (Clam) NEW NAME is LP - Limecola petalum 


#Tasks ####
#import the 2018 data and change the name of the clams to the new names




#Importing the data (Original Data Name = DE2018_alldata.csv)
meta18 <- read.csv("Oyster_data_raw/DE2018_alldata.csv")
meta18


#Changing the name to the new name
meta18$Species2 <- ifelse(meta18$Species =="CV", "CRASSOSTREA_VIRGINICA",
                            ifelse(meta18$Species == "IR", "ISCHADIUM_RECURVUM", 
                                   ifelse(meta18$Species == "MM", "AMERITELLA_MITCHELLI", "LIMECOLA_PETALUM")))


#Changing the abbreviations to the new name abbreviation
meta18$Species_abv_new <- ifelse(meta18$Species =="CV", "CV",
                          ifelse(meta18$Species == "IR", "IR", 
                                 ifelse(meta18$Species == "MM", "AM", "LP")))




#Oyster 16S Metadata Cleaning 
#Author: Diana Portugal 
#Contact: dportugal8@gmail.com 

#KEY ####
#Renaming the species collected
#CV - Crassostrea virginica (Oyster)
#IR - Ischadium recurvum (‘hooked’ Mussel)
#MM - Macoma mitchelli (Clam) NEW NAME is AM - Ameritella mitchelli
#MB = Macoma balthica (Clam) NEW NAME is LP - Limecola petalum 

library(tidyverse)
library(data.table)

#Importing the data (Original Data Name = metadata_de18.csv)
meta18 <- read.csv("Oyster_data_raw/metadata_de18.csv")
meta18

#Importing the data (Original Data Name = DE2018_alldata.csv)
de_data18 <- read.csv("Oyster_data_raw/DE2018_alldata.csv")
de_data18

#Importing the data (Original Data Name = asvtable_de18.csv)
asv18 <- fread("Oyster_data_raw/asvtable_de18.csv")


#Changing the name to the new name (DO NOT USE THIS FOR NOW)
#meta18$Species2 <- ifelse(meta18$Species =="CV", "CRASSOSTREA_VIRGINICA",
#ifelse(meta18$Species == "IR", "ISCHADIUM_RECURVUM", 
#ifelse(meta18$Species == "MM", "AMERITELLA_MITCHELLI", "LIMECOLA_PETALUM")))


#Changing the abbreviations to the new name abbreviation (META18)
meta18$Species <- ifelse(meta18$Species =="CV", "CV",
                          ifelse(meta18$Species == "IR", "IR", 
                                 ifelse(meta18$Species == "MM", "AM", "LP")))


#Changing the abbreviations to the new name abbreviation (DE_DATA18)
de_data18$Species <- ifelse(de_data18$Species =="CV", "CV",
                                 ifelse(de_data18$Species == "IR", "IR", 
                                        ifelse(de_data18$Species == "MM", "AM", "LP")))


#Creating a Treatment2 column with the combination ####
de_data18$Treatment2_18 <- paste(de_data18$Treatment1_Density, "_", de_data18$Treatment2_Diversity)

de_data18$Treatment2_18 <- ifelse(de_data18$Bucket =="HP1", "HIGH_POLY",
                                  ifelse(de_data18$Bucket =="HP2", "HIGH_POLY",
                                         ifelse(de_data18$Bucket =="HP3", "HIGH_POLY",
                                               ifelse(de_data18$Bucket =="HP4", "HIGH_POLY",
                                                     ifelse(de_data18$Bucket =="HP5", "HIGH_POLY",
                                                           ifelse(de_data18$Bucket =="HP6", "HIGH_POLY",
                                                                 ifelse(de_data18$Bucket =="HP7", "HIGH_POLY",
                                                                       ifelse(de_data18$Bucket =="HP8", "HIGH_POLY",
                                                                             ifelse(de_data18$Bucket =="HP9", "HIGH_POLY",
                                  
                            ifelse(de_data18$Bucket == "HM1", "HIGH_MONO",
                                   ifelse(de_data18$Bucket == "HM2", "HIGH_MONO",
                                          ifelse(de_data18$Bucket == "HM3", "HIGH_MONO",
                                                 ifelse(de_data18$Bucket == "HM4", "HIGH_MONO",
                                                        ifelse(de_data18$Bucket == "HM5", "HIGH_MONO",
                                                               ifelse(de_data18$Bucket == "HM6", "HIGH_MONO",
                                                                      ifelse(de_data18$Bucket == "HM7", "HIGH_MONO",
                                                                             ifelse(de_data18$Bucket == "HM8", "HIGH_MONO",
                                                                                    ifelse(de_data18$Bucket == "HM9", "HIGH_MONO",
                                   
                                   
                            ifelse(de_data18$Bucket == "LP1", "LOW_POLY",
                                   ifelse(de_data18$Bucket == "LP2", "LOW_POLY",
                                          ifelse(de_data18$Bucket == "LP3", "LOW_POLY",
                                                 ifelse(de_data18$Bucket == "LP4", "LOW_POLY",
                                                        ifelse(de_data18$Bucket == "LP5", "LOW_POLY",
                                                               ifelse(de_data18$Bucket == "LP6", "LOW_POLY",
                                                                      ifelse(de_data18$Bucket == "LP7", "LOW_POLY",
                                                                             ifelse(de_data18$Bucket == "LP8", "LOW_POLY",
                                                                                    ifelse(de_data18$Bucket == "LP9", "LOW_POLY", "LOW_MONO")))))))))))))))))))))))))))






#Creating the Column to merge the two data sets called merge18 & merge18
de_data18$Merge18 <- paste("2018",de_data18$Bucket, de_data18$Color.Number, de_data18$Species, sep="_") #(DE_DATA18)

meta18$Merge18 <- paste("2018", meta18$Color_Bucket, meta18$Number, meta18$Species, sep="_") #(META18)



#Using Merge to combine the two data frames  ####
#Done by matching the Merge18 columns
meta18data <- merge(meta18, de_data18, by = "Merge18", all.x = TRUE) #Matching by column Merge18, all.x referrers to Meta18 because it was on the X place



#Deleting columns in the new data frame
data_meta18_clean <- select(meta18data, 
                            -"X.x",
                            -"V1", 
                            -"Site",
                            -"peacrabs", 
                            -"Phase_1_DO", 
                            -"Phase_1_temp",
                            -"Phase_2_DO", 
                            -"Phase_2_Temp",
                            -"Overall_treatment", 
                            -"dead_barnacles",
                            -"Parasites", 
                            -"X.y", 
                            -"RFTM_score.y")


data_meta18_clean$delta_weight18 <- ((data_meta18_clean$Weight_post - data_meta18_clean$Weight)/data_meta18_clean$Weight)

#Creating a Treatment2 column with the combination ####
data_meta18_clean$Treatment2_18 <- paste(data_meta18_clean$Treatment1_Density, "_", data_meta18_clean$Treatment2_Diversity)



#Making Unique IDs the new row names for Phyloseq

write.csv(data_meta18_clean, file = "Oyster_data_raw/cleanmetadata18")

data_meta18_clean <- read.csv("Oyster_data_raw/cleanmetadata18")

rownames(data_meta18_clean) = data_meta18_clean$Merge18 

data_meta18_clean$Merge18=NULL

data_meta18_clean

write.csv(data_meta18_clean, file = "Oyster_data_raw/meta18cleaned")


#End here with the data cleaning and start a new script for the data analysis on phyloseq ####




























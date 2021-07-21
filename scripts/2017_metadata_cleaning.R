#Document Information ####
#Oyster 16S Metadata Cleaning 
#Author: Diana Portugal 
#Contact: dportugal8@gmail.com 


#Loading Data ####

#library(tidyverse)
#install.packages("data.table")
#library(data.table)

#Loading the data (Original Data is called DE_DATA_ForGenetics_17.csv)
de_data17 <- read.csv("Oyster_data_raw/DE_DATA_ForGenetics_17.csv")
de_data17

#Loading the data (Original Data is called metadata_de17.csv)
meta17 <- read.csv("Oyster_data_raw/metadata_de17.csv")
meta17

#Loading the data (Original Data Name = asvtable_de17.csv)
asv17 <- fread("Oyster_data_raw/asvtable_de17.csv")


#Renaming the Treatment Names ####
de_data17$Treatment
de_data17$Treatment2 <- ifelse(de_data17$Treatment =="HH", "HIGH_POLY",
                            ifelse(de_data17$Treatment == "HL", "HIGH_MONO", 
                                   ifelse(de_data17$Treatment == "LL", "LOW_MONO", "LOW_POLY")))


#Creating a new column names Colornumber 
de_data17$Colornumber <- paste0(de_data17$Color, de_data17$Number) 


#Creating the UniqueIDs in de_data ####
de_data17$UniqueID <- paste("2017", de_data17$Site, de_data17$Treatment2, de_data17$Colornumber, de_data17$Species, sep = "_")


#Using Merge to combine the two data frames  ####
#Done by matching the Unique ID columns
meta17data <- merge(meta17, de_data17, by = "UniqueID", all.x = TRUE) #Matching by column UniqueID, all.x referrers to Meta17 because it was on the X place


#Deleting columns in the new data frame
data_meta17_clean <- select(meta17data, 
                            -"X",
                            -"V1", 
                            -"Phase_1_DO",
                            -"Phase_1_temp", 
                            -"Phase_2_DO", 
                            -"Phase_2_Temp",
                            -"Overall_treatment", 
                            -"Date_post",
                            -"Notes_pre", 
                            -"POST_DEAD_ALIVE",
                            -"Dry_Weight_plate", 
                            -"Dry_weight_final", 
                            -"Dry_weight_shell", 
                            -"Notes_post", 
                            -"Genetics_Weight")




#Getting rid of the "missing" data ####
#Changing the MISSING data in the measurements to NA values

data_meta17_clean$Length_pre <-sub("MISSING","NA", data_meta17_clean$Length_pre)
data_meta17_clean

data_meta17_clean$Width_pre <-sub("MISSING","NA", data_meta17_clean$Width_pre)
data_meta17_clean

data_meta17_clean$Height_pre <-sub("MISSING","NA", data_meta17_clean$Height_pre)
data_meta17_clean

data_meta17_clean$Weight_pre <-sub("MISSING","NA", data_meta17_clean$Weight_pre)
data_meta17_clean



#Making Unique IDs the new row names for Phyloseq

write.csv(data_meta17_clean, file = "Oyster_data_raw/cleanmetadata17")

data_meta17_clean <- read.csv("Oyster_data_raw/cleanmetadata17")

rownames(data_meta17_clean) = data_meta17_clean$UniqueID

data_meta17_clean$UniqueID=NULL

data_meta17_clean

write.csv(data_meta17_clean, file = "Oyster_data_raw/meta17cleaned")



#End here with the data cleaning and start a new script for the data analysis on phyloseq ####


#DO NOT CHANGE ANYTHING HERE 









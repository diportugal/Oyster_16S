#Oyster 16S Metadata Cleaning 
#Author: Diana Portugal 
#Contact: dportugal8@gmail.com 

#KEY ####
#Renaming the species collected
#CV - Crassostrea virginica (Oyster)
#IR - Ischadium recurvum (‘hooked’ Mussel)
#MM - Macoma mitchelli (Clam) NEW NAME is AM - Ameritella mitchelli
#MB = Macoma balthica (Clam) NEW NAME is LP - Limecola petalum 


#Importing the data (Original Data Name = metadata_de18.csv)
meta18 <- read.csv("Oyster_data_raw/metadata_de18.csv")
meta18

#Importing the data (Original Data Name = DE2018_alldata.csv)
de_data18 <- read.csv("Oyster_data_raw/DE2018_alldata.csv")
de_data18

#Importing the data (Original Data Name = asvtable_de18.csv)
asv18 <- read.csv("Oyster_data_raw/asvtable_de18.csv")


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











filter(de_data18$Color_Bucket, is.na())


#Creating a new column called Bucket_Number
de_data18$Bucket_Number <- paste0(de_data18$Color_Bucket, de_data18$Number, is.na())

de_data17$Colornumber <- paste0(de_data17$Color, de_data17$Number) 


#Creating the UniqueID in DE_DATA18

de_data18$UniqueID <- paste("2018", de_data18$Treatment1_Density, de_data18$Treatment2_Diversity,  )











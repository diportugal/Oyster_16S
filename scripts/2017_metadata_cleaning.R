#Oyster 16S Metadata Cleaning 
#Author: Diana Portugal 
#Contact: dportugal8@gmail.com 

#This data set has the samples that were sequenced


#Loading the data (Original Data is called DE_DATA_ForGenetics_17.csv)
meta17 <- read.csv("Oyster_data_raw/DE_DATA_ForGenetics_17.csv")
meta17


#Renaming the Treatment Names 
meta17$Treatment
meta17$Treatment2 <- ifelse(meta17$Treatment =="HH", "HIGH_POLY",
                            ifelse(meta17$Treatment == "HL", "HIGH_MONO", 
                                   ifelse(meta17$Treatment == "LL", "LOW_MONO", "LOW_POLY")))

                       
                         
                          
                      
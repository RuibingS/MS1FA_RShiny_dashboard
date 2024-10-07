library(xcms)
library(ISFrag)
library(here)

# Feature Table Input - Si11 feature table
Si11_ft_xcms<-read.csv(here::here("Data","feature_table","Si11_feature_table_xcms.csv"))
colnames(Si11_ft_xcms)
Si11_ft_ISFrag<-Si11_ft_xcms[,-1]
colnames(Si11_ft_ISFrag)[1]<-"mz"
colnames(Si11_ft_ISFrag)[4]<-"rt"
head(Si11_ft_ISFrag)



# MS2 Annotation


# MS1directory specifies the full directory of the folder containing DDA all mzXML file(s).
MS2directory <- here::here("parameters_script","Si11_mzXML_files")


# The ms2.tofeaturetable() function assigns MS2 spectra from the provided DDA files to the MS1 feature table. It returns a new feature table with additional columns containing MS2 fragment information.
# Using XCMS feature table
featureTable <- ms2.assignment(MS2directory = MS2directory, XCMSFT = Si11_ft_ISFrag)


head(featureTable)

# Metabolite identification -> not working
# MOna library too big -> time out # directory containing the library file
# lib_name <- "MoNA-export-LC-MS-MS_Positive_Mode.msp" # name of the library file
# Si11.msp  <.not working <- use otehr function of read .msp ,if missing , NA 
# featureTable <- feature.annotation(featureTable = featureTable, lib_directory = lib_directory, lib_name = lib_name, dp = 0.8)
# head(featureTable)



##############################################################################################
##############################################################################################

# Part 4: Identification of ISF Features

# Identify level 3 in-source fragments.
# system.time(.level3 <- find.level3(MS1directory =here::here("parameters_script","Si11_mzXML_files"),
#                       MS1.files=list.files(here::here("parameters_script","Si11_mzXML_files")), 
#                       featureTable = featureTable, type = "multi"))
system.time(level3 <- find.level3(MS1directory =here::here("parameters_script","Si11_mzXML_files"),
                      MS1.files=here::here("parameters_script","Si11_mzXML_files","hm250 Si11mix 20uM C18 pos_GB5_01_35845.d.mzXML"), 
                      featureTable = featureTable, type = "single"))

str(level3)
level3[[10]]

# Identify level 2 in-source fragments.
level2 <- find.level2(ISFtable = level3)

# Identify level 1 in-source fragments.
level1 <- find.level1(ISF_putative = level2)


# Summarize ISFrag results after identifying all level 3, 2, and 1 in-source fragments.
results <- get.ISFrag.results(ISF_List = level1, featureTable = featureTable)

# Get complete feature table with all features and ISF relationship annotations.
resultFT <- export.ISFrag.results(ISFresult = results)
head(resultFT)


write.csv(resultFT,file ="~/Downloads/compare table/Si11_xcms_ISFrag_anno_V2.csv")




#############################################################################################################
#############################################################################################################
# ISFrag for antibiotics data set - PA14 feature table

# Feature Table Input

PA14_ft_xcms<-read.csv("~/Downloads/MS1FA_RShiny_dashboard/Data/feature_table/PA14_featureTable.csv")
colnames(PA14_ft_xcms)

# MS2 Annotation


# MS1directory specifies the full directory of the folder containing DDA mzXML file(s).
MS2directory <- "~/Downloads/MS1FA_RShiny_dashboard/Data/MS2_files/MS2_antibiotic"

# The ms2.tofeaturetable() function assigns MS2 spectra from the provided DDA files to the MS1 feature table. It returns a new feature table with additional columns containing MS2 fragment information.
# Using XCMS feature table
featureTable <- ms2.assignment(MS2directory = MS2directory, XCMSFT = PA14_ft_xcms,rt.tol = 20)

head(featureTable)


# Metabolite identification -> not working
# MOna library too big -> time out
# lib_directory <- "~/Downloads/MoNA" # directory containing the library file
# lib_name <- "MoNA-export-LC-MS-MS_Positive_Mode.msp" # name of the library file
# featureTable <- feature.annotation(featureTable = featureTable, lib_directory = lib_directory, lib_name = lib_name, dp = 0.8)
# head(featureTable)



##############################################################################################
##############################################################################################

# Part 4: Identification of ISF Features
MS1_files<-list.files("~/Downloads/PA14_mzXML_files/MS1_files_mzXML")
# Identify level 3 in-source fragments.
level3 <- find.level3(MS1directory = "~/Downloads/PA14_mzXML_files/MS1_files_mzXML",
                      MS1.files= MS1_files,
                      featureTable = featureTable, type = "multi")

# print - checking
str(level3)
level3[[10]]

# Identify level 2 in-source fragments.
level2 <- find.level2(ISFtable = level3)

# Identify level 1 in-source fragments.
level1 <- find.level1(ISF_putative = level2)


# Summarize ISFrag results after identifying all level 3, 2, and 1 in-source fragments.
results <- get.ISFrag.results(ISF_List = level1, featureTable = featureTable)

# Get complete feature table with all features and ISF relationship annotations.
resultFT <- export.ISFrag.results(ISFresult = results)
head(resultFT)


write.csv(resultFT,"~/Downloads/compare table/PA14_Antibiotics_xcms_ISFrag_V2.csv")































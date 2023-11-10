
rm(list=ls())

## set working directory
dir0 = "/Volumes/Bertolaz/HALO/funzioni_spaziale/Git_hub_scripts/"
dir0 = "/Volumes/Bertolaz/HALO/paper_Spaziale/Rmarkdown/"
setwd(dir0)

require(readxl)

### CLUSTER ANALYSIS
## we want to evaluate the cluster behavior of ACE2 cells
# null hypothesis: spatial random distribution
# alternative hypothesis: clustering behavior 

# data loading
Data = data.frame(read_xlsx("Supplementary_data.xlsx", sheet=1 ) )

## rename the variable names
names(Data)[which(names(Data)=="ACE2.Positive.classification")] ="ACE2"

## clustering test function
source("function_spatial_randomness_test.R")

res = spatial_randomness_test(data=Data,X = c("XMin", "XMax"), Y= c("YMin", "YMax"),  markers= c("ACE2"),
                          B=500, seed=123, scatterplot=TRUE,  colors=c("brown") )


######### SEGREGATION ANALYSIS
## we want to evaluate if AID and CD3 cells are significantly distant each other
# null hypothesis: spatial random distribution
# alternative hypothesis: segregation 

# data loading
Data = data.frame(read_xlsx("Supplementary_data.xlsx", sheet=2 ) )

## rename the variable names
names(Data)[which(names(Data)=="AID.Positive.Classification")] ="AID"
names(Data)[which(names(Data)=="CD3.Positive.Classification")] ="CD3"

## load function for segregation/aggregation analysis and running
source("function_spatial_dependence_test.R")

res =  spatial_dependence_test(data=Data,X=c("XMin", "XMax"), Y= c("YMin", "YMax"),  markers= c("AID", "CD3"),
                         B=5000, alternative = "segregation", colors = c( "red", "green"),  seed=1234)


## AGGREGATION ANALYSIS
## we want to evaluate if AID and GH2AX cells are significantly close each other
# null hypothesis: spatial random distribution
# alternative hypothesis: aggregation 
Data = data.frame(read_xlsx("Supplementary_data.xlsx", sheet=3 ) )

## rename the variable names
names(Data)[which(names(Data)=="aid.Positive.Classification")] ="AID"
names(Data)[which(names(Data)=="GH2AX.Positive.Nucleus.Classification")] ="GH2AX"

## we can consider the double positive cells as the third population.
# in this way we can study the spatial distribution of the double positive cell distribute 
# with respect to the single positive cells
Data$double = ifelse(Data$AID==1 & Data$GH2AX==1, 1, 0 )
Data$AID_only = ifelse(Data$AID==1 & Data$GH2AX==0, 1, 0)
Data$GH2X_only = ifelse(Data$AID==0 & Data$GH2AX==1, 1, 0)

## segregation/aggregation analysis function
source("function_spatial_dependence_test.R")

res =  spatial_dependence_test( data= Data,X= c("XMin", "XMax"), Y = c("YMin", "YMax"),  markers= c("AID_only", "GH2X_only", "double"),
                          B=500, alternative = "aggregation", colors = c("forestgreen", "orange", "red"),  seed=1234)


#############
######## EQUAL  DISTRIBUTION TEST 
Data = data.frame(read_xlsx("Supplementary_data.xlsx", sheet=4 ) )

# rename the variable names
names(Data)[which(names(Data)=="AID.Positive.Classification")] ="AID"
names(Data)[which(names(Data)=="CD3.Positive.Classification")] ="CD3"

# loading the function for distribution equality testing
source("function_distribution_equality_test.R")

res = distribution_equality_test(data=Data, X= c("XMin", "XMax") , Y= c("YMin", "YMax"),  
                                 M1="CD3" , M2="AID", colM1="green" ,  colM2="red" ,
                                 rm.double=TRUE, B=5000)


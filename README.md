# SpaceR
Spatial Cell Resampling for quantitative analysis of cytometry data

<style type="text/css">
.main-container {
  max-width: 1800px !important;
  margin-left: auto;
  margin-right: auto;
}
</style>

![](README_files/figure-markdown_github/Cellula_spaziale.png)

# Overview

The **SpaceR** package implement point pattern spatial analysis for
spatial cytometry data.

The SpaceR functions allow running novel resampling methods for testing
spatial randomness and evaluating relationships among different cell
populations. Our methods obviate the need for domain area estimation and
provide non-biased clustering measurements.

In this tutorial, the input file contains the tissue cell coordinates
and cell markers information. We will evaluate the presence of cell
clustering behavior, the presence of aggregation and segregation effects
between cell populations, and the differences in spatial cell
distributions.

Our approach can be applied to any point pattern analysis in which the
study of sub-populations is limited by an inaccurate estimation of the
domain area.

**Reference**: Bertolazzi, Tumminello, Morello, Belmonte, Tripodo (2023)
*Resampling approaches for the quantitative analysis of spatially
distributed cells*

# Getting Started

**Downloadable files**:

The input files and functions to run the following analyzes can be
downloaded from the SpaceR GitHub page:

-   *Supplementary_data.xlsx*: input data used in this tutorial
-   *functions_description.pdf*: detailed description of function
    parameters
-   *function_spatial_randomness_test.R*: function which performs the
    spatial randomness test
-   *function_spatial_dependence_test.R*: function which performs the
    spatial dependence
-   *function_distribution_equality_test.R*: function which performs
    distribution equality test
-   *run_functions.R*: R script which show how to run the functions

**Input data**:

-   *X-axis and Y-axis cell coordinates*
-   *Marker positivities that identify cell populations*

**Functions to import into R enviroment**:

-   *spatial_randomness_test*: function which performs the spatial
    randomness test
-   *spatial_dependence_test*: function which performs the spatial
    dependence
-   *distribution_equality_test*: function which performs distribution
    equality test

**R packages**:

The list of packages we are going to use for this analysis:

``` r
  library(readxl)
  library(usedist)
  library(TSdist)
  library(graphics)
  require(spatstat)
```

# Usage

## Spatial randomness test

The spatial randomness test evaluates the presence of cell clustering
behavior. In this example we study the possible cluster behavior of the
ACE2 cell population.

Data loading:

``` r
Data = data.frame(read_xlsx("Supplementary_data.xlsx", sheet=1))

head(Data)
```

    ##    Image.File.Name Analysis.Region                     Algorithm.Name Object.Id
    ## 1 Snap-111803.tiff    entire image Indica Labs - Multiplex IHC v3.0.3         0
    ## 2 Snap-111803.tiff    entire image Indica Labs - Multiplex IHC v3.0.3         1
    ## 3 Snap-111803.tiff    entire image Indica Labs - Multiplex IHC v3.0.3         2
    ## 4 Snap-111803.tiff    entire image Indica Labs - Multiplex IHC v3.0.3         3
    ## 5 Snap-111803.tiff    entire image Indica Labs - Multiplex IHC v3.0.3         4
    ## 6 Snap-111803.tiff    entire image Indica Labs - Multiplex IHC v3.0.3         5
    ##   XMin XMax YMin YMax ACE2.Positive.classification ACE2.Nucleus.OD
    ## 1 1862 1892  797  840                            1        0.348051
    ## 2 1889 1916  798  823                            1        0.411348
    ## 3 1828 1851  799  833                            0        0.037044
    ## 4 1682 1707  801  829                            1        0.275323
    ## 5 1788 1817  802  837                            0        0.103893
    ## 6 1764 1790  806  834                            1        0.248718
    ##   Cell.Area..µm.. Cytoplasm.Area..µm.. Nucleus.Area..µm..
    ## 1          524768               280899             243869
    ## 2          267674               184092              83582
    ## 3           33327               220064             113206
    ## 4           29095                19573               9522
    ## 5          345966               235934             110032
    ## 6          305762               200491             105271
    ##   Nucleus.Perimeter..µm. Nucleus.Roundness
    ## 1                  24.84          0.621809
    ## 2                  12.65          0.848912
    ## 3                  15.18          0.684204
    ## 4                  13.34          0.814121
    ## 5                   18.4          0.563253
    ## 6                   13.8          0.839799

This dataset contains cytometric information such as cell area and
perimeter. To run the spatial analysis we only need cell coordinates and
marker classifications.

The *XMin* and *XMax* variables contain the X-axis coordinates of cell
perimeters. If we provide both variables as X-axis input, the function
for testing the spatial randomness automatically calculates the average
of X-axis coordinates. This average approximately corresponds to the
X-axis cell centroid. Same type of input is considered for the Y-axis
coordinates.

The *ACE2.Positive.classification* variable contains the cell marker
classification. For simplicity, we rename this variable:

``` r
names(Data)[which(names(Data)=="ACE2.Positive.classification")]="ACE2"                               
```

Loading and running the function to test the spatial randomness:

``` r
source("function_spatial_randomness_test.R")

res = spatial_randomness_test(data=Data,
                          X = c("XMin", "XMax"), 
                          Y= c("YMin", "YMax"),
                          markers= c("ACE2"),
                          B=5000, seed=123, 
                          scatterplot=TRUE, 
                          colors=c("brown") )
```

    ## Spatial randomness testing... 
    ## 20%  Done  
    ## 40%  Done  
    ## 60%  Done  
    ## 80%  Done  
    ## 100%  Done  
    ##  
    ##        n average_NN rand_average_NN     R   pval    
    ## ACE2 203      47.34           61.35 0.772 <2e-04 ***
    ## --- 
    ## alternative hypothesis: spatial clustering
    ## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

![](README_files/figure-markdown_github/SpatialRandomness-1.png)

Following are the function arguments:

-   data: dataframe containing the cell spatial coordinates and marker
    positivities (1 = positive, 0 = negative). Negative cells must be
    included in the dataset
-   X: vector containing the names of variables that indicate the
    spatial coordinates of cells along the X-axis. The average values
    will be considered in the case of two X-axis variables
-   Y: vector containing the names of variables that indicate the
    spatial coordinates of cells along the Y-axis. The average values
    will be considered in the case of two Y-axis variables
-   markers: vector containing the names of variables that indicate
    marker positivities
-   B: number of random sample extractions
-   seed: seed of the resampling
-   scatterplot: logicals. If TRUE, the scatterplot of the cell spatial
    distributions will be plotted
-   colors: vector of colours for the scatterplot representation. The
    colour order should be the same as the marker vector

Function output:

``` r
res
```

    ##   marker   n average_NN rand_average_NN     R   pval
    ## 1   ACE2 203      47.34           61.35 0.772 <2e-04

where,

-   marker: cell population under exam
-   n: number of cells which compose the cell population
-   average_NN: average nearest neighbor distance inside the cell
    population
-   rand_average_NN: average nearest neighbor distance inside the random
    populations of size *n*
-   R: empirical *R* index. *R* \< 1 indicates aggregation (i.e.,
    clustering), *R* \> 1 indicates segregation, *R* ≈ 1 indicates
    spatial randomness
-   pval: spatial randomness *p*-value. Null hypothesis: spatial
    randomness. Alternative hypothesis: cell aggregation

The significant results we observed indicate the presence of spatial
clustering in the ACE2 cell population.

## Segregation analysis

The spatial dependence test can be used to evaluate the presence of
segregation or aggregation between two cell populations. In this
example, we evaluate if AID and CD3 cells significantly segregate.

``` r
Data = data.frame(read_xlsx("Supplementary_data.xlsx", sheet=2 ) )

names(Data)[which(names(Data)=="AID.Positive.Classification")] ="AID"
names(Data)[which(names(Data)=="CD3.Positive.Classification")] ="CD3"

source("function_spatial_dependence_test.R")

res =  spatial_dependence_test(data=Data, 
                               X=c("XMin", "XMax"), 
                               Y= c("YMin", "YMax"), 
                               markers= c("AID", "CD3"),
                               B=5000, seed=1234, 
                               alternative = "segregation", 
                               colors = c( "red", "green"))
```

    ## Cell resampling... 
    ## 20%  Done  
    ## 40%  Done  
    ## 60%  Done  
    ## 80%  Done  
    ## 100%  Done  
    ## --- 
    ## Average distances of nearest neighbor cells to AID (n =  199 ):  
    ##       n average_NN rand_average_NN     R   pval    
    ## CD3 158     265.84          100.48 2.646 <2e-04 ***
    ## --- 
    ## alternative hypothesis: segregation 
    ## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

![](README_files/figure-markdown_github/segregation-1.png)

In the *marker* vector given as fuction input, AID is the baseline
population from which the nearest neighbor distances are calculated. In
this example, for each AID cell we search the nearest CD3 cell and
calculate the distance.

``` r
res
```

    ##   markers   n average_NN rand_average_NN     R   pval
    ## 1     CD3 158     265.84          100.48 2.646 <2e-04

Output:

-   average_NN: average distance between AID and CD3 populations
-   rand_average_NN: average distance between AID cells and the random
    populations
-   R: empirical *R*<sub>*a**b*</sub> index. *R* \< 1 indicates
    aggregation effect between the two populations, *R* \> 1 indicates
    segregation effect between the two populations, *R* ≈ 1 indicates
    the absence of spatial dependence between the two populations
-   pval: spatial dependence *p*-value. Null hypothesis: spatial
    independence. Alternative hypothesis: cell segregation

## Aggregation analysis

The spatial dependence test can also be used to evaluate the presence of
aggregation between two cell populations. In this example, we evaluate
if AID and GH2AX cells are significantly close to each other.

We have observed double positive cells. A possible analysis strategy is
to consider them as an additional population. In this way we can assess
how double positive cells distribute and interact with the other cell
populations.

``` r
Data = data.frame(read_xlsx("Supplementary_data.xlsx", sheet=3 ) )

names(Data)[which(names(Data)=="aid.Positive.Classification")] ="AID"
names(Data)[which(names(Data)=="GH2AX.Positive.Nucleus.Classification")] ="GH2AX"

Data$double = ifelse(Data$AID==1 & Data$GH2AX==1, 1, 0 )
Data$AID_only = ifelse(Data$AID==1 & Data$GH2AX==0, 1, 0)
Data$GH2X_only = ifelse(Data$AID==0 & Data$GH2AX==1, 1, 0)

source("function_spatial_dependence_test.R")

res =  spatial_dependence_test( data= Data,
                                X= c("XMin", "XMax"), 
                                Y = c("YMin", "YMax"),  
                                markers= c("AID_only", "GH2X_only", "double"),
                                B=5000,   seed=1234,
                                alternative = "aggregation", 
                                colors = c("forestgreen", "orange", "red") )
```

    ## Cell resampling... 
    ## 20%  Done  
    ## 40%  Done  
    ## 60%  Done  
    ## 80%  Done  
    ## 100%  Done  
    ## --- 
    ## Average distances of nearest neighbor cells to AID_only (n =  988 ):  
    ##             n average_NN rand_average_NN     R   pval    
    ## GH2X_only 163      72.14           91.61 0.787 <2e-04 ***
    ## double    102     100.03          112.63 0.888 0.0936   .
    ## --- 
    ## alternative hypothesis: aggregation 
    ## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

![](README_files/figure-markdown_github/aggregation-1.png)

``` r
res
```

    ##     markers   n average_NN rand_average_NN     R   pval
    ## 1 GH2X_only 163      72.14           91.61 0.787 <2e-04
    ## 2    double 102     100.03          112.63 0.888 0.0936

As we shown in the previous example, AID is the baseline population from
which the nearest neighbor distances are calculated (AID was given as
the first element of the *marker* vector). Therefore, each row of the
output matrix describes the relation between AID cells and an other cell
population. The results suggest the presence of aggregation behavior
between AID and GH2X cells, while the double positive cells don’t
significantly aggregate with AID cells.

## Distribution equality test

The distribution equality test evaluates if two cell populations are
differentially distributed over the domain. In this example, we compare
AID and CD3 cell populations.

``` r
Data = data.frame(read_xlsx("Supplementary_data.xlsx", sheet=4 ) )

names(Data)[which(names(Data)=="AID.Positive.Classification")] ="AID"
names(Data)[which(names(Data)=="CD3.Positive.Classification")] ="CD3"

source("function_distribution_equality_test.R")

res = distribution_equality_test(data=Data, 
                                 X= c("XMin", "XMax"),
                                 Y= c("YMin", "YMax"),  
                                 M1="CD3" , M2="AID", 
                                 colM1="green" ,  colM2="red",
                                 rm.double=TRUE, 
                                 B=5000)
```

    ## warning:  2  double positive cells have been removed 
    ##  
    ##  Distribution equality testing... 
    ## 20%  Done  
    ## 40%  Done  
    ## 60%  Done  
    ## 80%  Done  
    ## 100%  Done  
    ## 

![](README_files/figure-markdown_github/equality-1.png)

    ## test.statistic =  1679.62 , p-value  < 2e-04  
    ## alternative hypothesis:  CD3  and  AID cells have different spatial distributions 
    ## dissimilarity index =  12.03

The results indicate a different cell spatial distribution between the
two populations. The left part of the graph shows the cell distribution,
while the right part compares the expected and the observed frequency
differences.

# Invasive buffelgrass, Cenchrus ciliaris, balances opportunistic acquisition of foliar fungi with host and environmental filtering in its introduced range.

---
github respository name: Bowmanetal2023_BuffelGrass_FoliarFungalCommunities
---
author: Dr. Liz Bowman
---
contact: eabowman@utexas.edu
---
Data and analysis code associated with the manuscript "Invasive buffelgrass, Cenchrus ciliaris, balances opportunistic acquisition of foliar fungi with host and environmental filtering in its introduced range."

The following code will load libraries and run the analyses. Output files will go into folders marked figures/ and results/. Please create these folders prior to running the code. Output names should match figure and table names in the manuscript.

## Load libraries
```{r, include = F}
source("LoadLibraries.R")
```

## Question 1: How does an invasive plant’s foliar fungal community overlaps between its native and introduced ranges?
```{r}
source("Q1_NativeIntroducedRange_CCiliaris.R")
```

## Question 2: What drives foliar fungal community composition within the introduced range compared to co-occurring native and non-native plants?
```{r}
source("Q2_IntroducedRange_ComparisonCoOccurring.R")
```

## Question 3: Are foliar fungi associated with C. ciliaris are co-introduced or locally assembled?
```{r}
source("Q3_OriginOfSymbiontsAssociatedwithCCiliaris.R")
```

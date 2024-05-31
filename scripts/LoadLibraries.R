# This script is to load all libraries required for analsyes in the following
# files:
# 1. Q1_NativeIntroducedRange_CCiliaris.R
# 2. Q2_IntroducedRange_ComparisonCoOccurring.R
# 3. Q3_OriginOfSymbiontsAssociatedwithCCiliaris.R
# 
# See the ReadMe.R file for an overview of scripts and data

# install.packages('tidyverse')
library(tidyverse)

# install.packages('vegan')
library(vegan)

# install.packages('ape')
library(ape)

# install.packages("indicspecies")
library(indicspecies)

# install.packages("cooccur")
library(cooccur)

# https://github.com/traitecoevo/fungaltraits/
# install.packages("devtools")
# devtools::install_github("ropenscilabs/datastorr")
# devtools::install_github("traitecoevo/fungaltraits")
library(fungaltraits)

install.packages("remotes")
remotes::install_github("kassambara/rstatix")

install.packages("rcompanion")

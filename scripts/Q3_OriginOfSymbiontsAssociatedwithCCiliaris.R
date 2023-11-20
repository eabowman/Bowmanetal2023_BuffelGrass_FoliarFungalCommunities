## Script created by Dr. Liz Bowman (eabowman@utexas.edu), May 14, 2021 
## Focus is to answer question 3: To what extent are foliar fungi associated
## with C. ciliaris are co-introduced or locally assembled.
## 
## Data: 95% Sequence Similarity Rarefied data (SitexSpecies_95sim_Raref.csv)
## There are two parts:
## Part 1 uses the C. cilaris only data from both the native and introduced range
## A: Isolate part 1 data
## B: Indicator species and co-occurrence patterns
## C: Fungal Trophic Mode
## 
## Part 2 uses data from all species in the introduced range only
## D: Isolate part 2 data
## E: Dissimilarity comparison by host
## F: Taxonomic groups shared by all grass species in introduced range
## G: Comparison of OTU shared between ranges with OTU found in introduced range

# Be sure to run the LoadLibraries.R file before running code below.

#========================================================================================#
# Part 1: C. cilaris data from native and introduced range ----------------
#========================================================================================#

# read in data
otu.data <- read.csv('E13-133.data.output/SitexSpecies_95sim_Raref.csv')
clim.data <- read_csv('E13-133.data.output/ClimateData.csv')
other.data <- read_csv('E13-133.data.output/fungal.taxonomy/eukaryote_zotu95-Phylum.csv')

#========================================================================================#
## A: Isolate data from C. ciliaris hosts in TX, AZ, and Kenya----
#========================================================================================#

# Removed data from other plant parts than leaves (roots, dead leaves, litter)
# Removed data from plants that were washed or placed in ethanol
otu.data %>%
  filter(Host %in% c("Cenchrus_ciliaris"),
         Plant_Part == "Leaf Green",
         Wash == "Unwashed") -> buffel.data

# remove non-fungal Otus
buffel.data <- buffel.data[!colnames(buffel.data) %in% other.data$Otu]

# isolate same data from clim.data
clim.data <- clim.data[clim.data$Illumina %in% buffel.data$samples, ]

### add metadata to buffel data ----
for(i in buffel.data$samples){
  # prec.i <- meta.data[meta.data$Illumina == i, 'Rain_mm_30yr_Av']
  # temp.i <- meta.data[meta.data$Illumina == i, 'Tmean_C_30yr']
  # buffel.data[buffel.data$samples == i, 'MeanPrecip'] <- prec.i
  # buffel.data[buffel.data$samples == i, 'MeanTemp'] <- temp.i
  prec.18.i <- clim.data[clim.data$Illumina == i, 'Prec.18']
  tmax.18.i <- clim.data[clim.data$Illumina == i, 'Tmax.18']
  ann.prec.hist.i <- clim.data[clim.data$Illumina == i, 'BIO12'] # Annual mean prec
  wet.prec.hist.i <- clim.data[clim.data$Illumina == i, 'BIO13'] # Prec wettest quarter
  warm.temp.hist.i <- clim.data[clim.data$Illumina == i, 'BIO10'] # mean Temp Warmest quarter
  ann.temp.hist.i <- clim.data[clim.data$Illumina == i, 'BIO1'] # Annual mean temp
  lat.i <- clim.data[clim.data$Illumina == i, 'lat']
  long.i <- clim.data[clim.data$Illumina == i, 'lat']
  state.i <- clim.data[clim.data$Illumina == i, 'Location.state']
  buffel.data[buffel.data$samples == i, 'Prec.18'] <- prec.18.i
  buffel.data[buffel.data$samples == i, 'Tmax.18'] <- tmax.18.i
  buffel.data[buffel.data$samples == i, 'MAP'] <- ann.prec.hist.i
  buffel.data[buffel.data$samples == i, 'MAT'] <- ann.temp.hist.i
  buffel.data[buffel.data$samples == i, 'MPWQ'] <- wet.prec.hist.i
  buffel.data[buffel.data$samples == i, 'MTWQ'] <- warm.temp.hist.i
  buffel.data[buffel.data$samples == i, 'lat'] <- lat.i
  buffel.data[buffel.data$samples == i, 'long'] <- long.i
  buffel.data[buffel.data$samples == i, 'Location.state'] <- state.i
}

# rearrange
buffel.data <- buffel.data[c(1:10, 357:365, 11:356)]

## PCNM----
# Geographic distance
geo.dist <- vegdist(buffel.data[c('lat','long')], method = 'euclidean')
geo.pcnm <- pcnm(geo.dist)
### add to buffel.data----
buffel.data['Geo.pcnm.1'] <- geo.pcnm$vectors[,1]
buffel.data['Geo.pcnm.2'] <- geo.pcnm$vectors[,2]
# reorder
buffel.data <- buffel.data[c(1:19, 366:367, 20:365)]

# Remove columns with no occurrences (OTUs unique to other species than buffel grass)
comm.data <- buffel.data[22:length(buffel.data)]
comm.data[colSums(comm.data) > 0] -> comm.data
cbind(buffel.data[1:21], comm.data) -> buffel.data

#========================================================================================#
## B: Indicator and Co-occurence analyses----
#========================================================================================#

# import taxonomic data
tax.data <- read.csv('E13-133.results/OTU_95SeqSim_Rarefied/Native.IntroducedRange.Shared.Unique.Otu/Overlap.tax.csv')

## Indicator species ----
# isolate community data
buffel.comm <- buffel.data[22:length(buffel.data)] 

# change table to presence absence (0s and 1s)
  for(r in 1:nrow(buffel.comm)){
    for(c in 1:ncol(buffel.comm)){
      if(buffel.comm[r,c] >= 1){
        buffel.comm[r,c] <- 1}
      else {buffel.comm[r,c] <- 0}}}

# remove columns with < 1 occurrences
buffel.comm[colSums(buffel.comm) > 1] -> buffel.comm

# indicator analysis assessing what species are indicative of range
ind.range <- multipatt(buffel.comm,
                       cluster = buffel.data$Range,
                       duleg = F,
                       func = 'r.g',
                       control = how(nperm = 999))

# filter Otu with p < 0.05
ind.range$sign[ind.range$sign$p.value < 0.05,] -> ind.range.sig

# filter out NA
filter(ind.range.sig, p.value > 0) -> ind.range.sig

# add taxonomic data
# First add Otus as a column instead of just row names
ind.range.sig$Otu <- row.names(ind.range.sig)
left_join(ind.range.sig, tax.data) -> tax.ind.range.sig

# write.csv(tax.ind.range.sig,
#           'Table3.csv',
#           row.names = F)


## Cooccur analyses ----
### Native range ----
# isolate Kenya data
filter(buffel.data, Range == 'Native') -> buffel.nat.data

# isolate community data
buffel.nat.data[22:length(buffel.nat.data)] -> buffel.nat.comm

# change table to 1s and 0s (presence/absence)
for(r in 1:nrow(buffel.nat.comm)){
  for(c in 1:ncol(buffel.nat.comm)){
    if(buffel.nat.comm[r,c] >= 1){
      buffel.nat.comm[r,c] <- 1}
    else {buffel.nat.comm[r,c] <- 0}
  }
}

# remove columns (OTUs) with no observations (colSums = 0) or only 1 observation
buffel.nat.comm[colSums(buffel.nat.comm) > 1] -> buffel.nat.comm

# Transpose the table
buffel.nat.comm.t <- t(buffel.nat.comm)

cooccur.nat.buffel <- cooccur(buffel.nat.comm.t,
                          type = 'spp_site',
                          thresh = T,
                          spp_names = T)

summary(cooccur.nat.buffel)
plot(cooccur.nat.buffel)

# Return a table of results
as.data.frame(prob.table(cooccur.nat.buffel)) -> prob.nat

# isolate co-occurrences that are sig. positive or negative.
prob.nat %>%
  filter(p_gt < 0.05) -> prob.nat.pos.sig

prob.nat %>%
  filter(p_lt < 0.05) -> prob.nat.neg.sig

rbind(prob.nat.pos.sig, prob.nat.neg.sig) -> prob.nat.sig

# add taxonomic data
prob.nat.sig$sp1_tax_class <- NA
prob.nat.sig$sp2_tax_class <- NA
prob.nat.sig$sp1_tax_family <- NA
prob.nat.sig$sp2_tax_family <- NA
prob.nat.sig$sp1_tax_genus <- NA
prob.nat.sig$sp2_tax_genus <- NA

for(i in unique(prob.nat.sig$sp1_name)){
  if(i %in% tax.data$Otu){
  tax.i <- tax.data[tax.data$Otu == i, 'Class']
  tax.f.i <- tax.data[tax.data$Otu == i, 'Family']
  tax.g.i <- tax.data[tax.data$Otu == i, 'Genus']
  prob.nat.sig[prob.nat.sig$sp1_name == i, 'sp1_tax_class'] <- tax.i
  prob.nat.sig[prob.nat.sig$sp1_name == i, 'sp1_tax_family'] <- tax.f.i
  prob.nat.sig[prob.nat.sig$sp1_name == i, 'sp1_tax_genus'] <- tax.g.i
  }
}

for(i in unique(prob.nat.sig$sp2_name)){
  if(i %in% tax.data$Otu){
  tax.i <- tax.data[tax.data$Otu == i, 'Class']
  tax.f.i <- tax.data[tax.data$Otu == i, 'Family']
  tax.g.i <- tax.data[tax.data$Otu == i, 'Genus']
  prob.nat.sig[prob.nat.sig$sp2_name == i, 'sp2_tax_class'] <- tax.i
  prob.nat.sig[prob.nat.sig$sp2_name == i, 'sp2_tax_family'] <- tax.f.i
  prob.nat.sig[prob.nat.sig$sp2_name == i, 'sp2_tax_genus'] <- tax.g.i
  }
}

# write.csv(prob.nat.sig,
#           'CoOccur_NativeRange.csv',
#           row.names = F)

prob.nat.sig %>%
  filter(p_gt < 0.05) %>%
  group_by(sp1_tax_class, sp2_tax_class) %>%
  summarize(count = n()) -> assess.nat

# plot of expected and observed co-occurrence as a function of positive, negative,
# and random cooccurrences. 
# add positive, negative, and random column to prob.nat 
for(i in 1:nrow(prob.nat)){
  if(prob.nat[i, 'p_gt' ] < 0.05){
    prob.nat[i, 'cooccurrence'] <- 'Positive'
  } else if(prob.nat[i, 'p_lt'] < 0.05) {
    prob.nat[i, 'cooccurrence'] <- 'Negative'
  } else {prob.nat[i, 'cooccurrence'] <- 'Random'}
}
main.nat <- ggplot(prob.nat, aes(x = exp_cooccur,
                     y = obs_cooccur,
                     fill = cooccurrence,
                     group = 1)) +
  geom_point(size = 2,
             shape = 21,
             color = 'black',
             alpha = 0.8,
             position = position_jitter(h = 0.3,
                                        w = 0.3)) +
  scale_fill_manual(values = c('red', 'blue', 'grey')) +
  geom_smooth(method = 'lm',
              se = F,
              formula = y ~ x,
              color = 'black',
              linetype = 'dotdash') +
  ylab('Observed co-occurrences') +
  xlab('Expected co-occurrences') +
  theme_classic() +
  theme(axis.text = element_text(size = 18, color = 'black'),
        axis.title = element_text(size = 20),
        legend.position = 'none')


# Assess proportion that is random versus non-random
prob.nat %>%
  group_by(cooccurrence) %>%
  summarise(count = n()) %>%
  mutate(percentage = count/sum(count)*100) -> summary.prob.nat

inset.nat <- ggplot(summary.prob.nat,
       aes(y = percentage,
           x = cooccurrence,
           fill = cooccurrence)) +
  geom_bar(stat = 'identity',
           color = 'black',
           alpha = 0.6) +
  xlab('') +
  ylab('') +
  scale_fill_manual(values = c('red', 'blue', 'grey')) +
  theme_classic() +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14),
        legend.position = 'none')

cowplot::ggdraw() +
  cowplot::draw_plot(main.nat) +
  cowplot::draw_plot(inset.nat, x = 0.55, y = 0.14,
                     width = 0.4, height = 0.4) -> full.nat

# ggsave('Fig5A.jpeg',
#        plot = full.nat, device = 'jpeg',
#        width = 7, height = 5, unit = 'in')


### Introduced range ----
#### All samples -----
# isolate US data
filter(buffel.data, Range != 'Native') -> buffel.int.data

# isolate community data
buffel.int.data[22:length(buffel.int.data)] -> buffel.int.comm

# change table to 1s and 0s (presence/absence)
for(r in 1:nrow(buffel.int.comm)){
  for(c in 1:ncol(buffel.int.comm)){
    if(buffel.int.comm[r,c] >= 1){
      buffel.int.comm[r,c] <- 1}
    else {buffel.int.comm[r,c] <- 0}
  }
}

# remove columns (OTUs) with no observations (colSums = 0) and with one observation
buffel.int.comm[colSums(buffel.int.comm) > 1] -> buffel.int.comm

buffel.int.comm.t <- t(buffel.int.comm)

cooccur.int.buffel <- cooccur(buffel.int.comm.t,
                          type = 'spp_site',
                          thresh = T,
                          spp_names = T)

summary(cooccur.int.buffel)
plot(cooccur.int.buffel)

# Return a table of results
as.data.frame(prob.table(cooccur.int.buffel)) -> prob.int

prob.int %>%
  filter(p_gt < 0.05) -> prob.int.pos.sig

prob.int %>%
  filter(p_lt < 0.05) -> prob.int.neg.sig

rbind(prob.int.pos.sig, prob.int.neg.sig) -> prob.int.sig

# add taxonomic data
prob.int.sig$sp1_tax_class <- NA
prob.int.sig$sp2_tax_class <- NA
prob.int.sig$sp1_tax_family <- NA
prob.int.sig$sp2_tax_family <- NA
prob.int.sig$sp1_tax_genus <- NA
prob.int.sig$sp2_tax_genus <- NA

for(i in unique(prob.int.sig$sp1_name)){
  if(i %in% tax.data$Otu){
  tax.i <- tax.data[tax.data$Otu == i, 'Class']
  tax.f.i <- tax.data[tax.data$Otu == i, 'Family']
  tax.g.i <- tax.data[tax.data$Otu == i, 'Genus']
  prob.int.sig[prob.int.sig$sp1_name == i, 'sp1_tax_class'] <- tax.i
  prob.int.sig[prob.int.sig$sp1_name == i, 'sp1_tax_family'] <- tax.f.i
  prob.int.sig[prob.int.sig$sp1_name == i, 'sp1_tax_genus'] <- tax.g.i
  }
}

for(i in unique(prob.int.sig$sp2_name)){
  if(i %in% tax.data$Otu){
  tax.i <- unique(tax.data[tax.data$Otu == i, 'Class'])
  tax.f.i <- tax.data[tax.data$Otu == i, 'Family']
  tax.g.i <- tax.data[tax.data$Otu == i, 'Genus']
  prob.int.sig[prob.int.sig$sp2_name == i, 'sp2_tax_class'] <- tax.i
  prob.int.sig[prob.int.sig$sp2_name == i, 'sp2_tax_family'] <- tax.f.i
  prob.int.sig[prob.int.sig$sp2_name == i, 'sp2_tax_genus'] <- tax.g.i
  }
}

# write.csv(prob.int.sig,
#           'CoOccur_IntroducedRange.csv',
#           row.names = F)

prob.int.sig %>%
  filter(p_gt < 0.05) %>%
  group_by(sp1_tax_class, sp2_tax_class) %>%
  summarize(count = n()) -> assess.int

# plot of expected and observed co-occurrence as a function of positive, negative,
# and random cooccurrences. 
# add positive, negative, and random column to prob.nat 
for(i in 1:nrow(prob.int)){
  if(prob.int[i, 'p_gt' ] < 0.05){
    prob.int[i, 'cooccurrence'] <- 'Positive'
  } else if(prob.int[i, 'p_lt'] < 0.05) {
    prob.int[i, 'cooccurrence'] <- 'Negative'
  } else {prob.int[i, 'cooccurrence'] <- 'Random'}
}

main.int <- ggplot(prob.int, aes(x = exp_cooccur,
                     y = obs_cooccur,
                     fill = cooccurrence,
                     group = 1)) +
  geom_point(size = 2,
             shape = 21,
             color = 'black',
             alpha = 0.6,
             position = position_jitter(h = 0.3,
                                        w = 0.3)) +
  scale_fill_manual(values = c('red', 'blue', 'grey')) +
  xlim(0,50) +
  geom_smooth(method = 'lm',
              se = F,
              formula = y ~ x,
              color = 'black',
              linetype = 'dotdash') +
  ylab('Observed co-occurrences') +
  xlab('Expected co-occurrences') +
  theme_classic() +
  theme(axis.text = element_text(size = 18, color = 'black'),
        axis.title = element_text(size = 20),
        legend.position = 'none')

# Assess proportion that is random versus non-random
prob.int %>%
  group_by(cooccurrence) %>%
  summarise(count = n()) %>%
  mutate(percentage = count/sum(count)) -> summary.prob.int

inset.int <- ggplot(summary.prob.int,
       aes(y = percentage,
           x = cooccurrence,
           fill = cooccurrence)) +
  geom_bar(stat = 'identity', color = 'black', alpha = 0.6) +
  xlab('') +
  ylab('') +
  scale_fill_manual(values = c('red', 'blue', 'grey')) +
  theme_classic() +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14),
        legend.position = 'none')

cowplot::ggdraw() +
  cowplot::draw_plot(main.int) +
  cowplot::draw_plot(inset.int, x = 0.55, y = 0.14,
                     width = 0.4, height = 0.4) -> full.int

# ggsave('Fig5B.jpeg',
#        plot = full.int, device = 'jpeg',
#        width = 7, height = 5, unit = 'in')


#### Randomly subsampled 27 samples to match number from native range -----
# isolate US data
filter(buffel.data, Range != 'Native') -> buffel.int.data

results.frame <- data.frame(rep = 1:100,
                            total.pairs = NA, 
                            positive = NA,
                            positive.percent = NA,
                            negative = NA,
                            negative.percent = NA)

random.introd <- function(dataframe, index.num){
  buffel.sub <- dataframe[sample(1:nrow(dataframe), 27, replace = F), ]
  # isolate community data
  buffel.comm <- buffel.sub[22:length(buffel.sub)]
  # change table to presence absence (0s and 1s)
  for(r in 1:nrow(buffel.comm)){
    for(c in 1:ncol(buffel.comm)){
      if(buffel.comm[r,c] >= 1){
        buffel.comm[r,c] <- 1}
      else {buffel.comm[r,c] <- 0}}}
  # remove columns (OTUs) with no observations (colSums = 0) and with one observation
  buffel.comm[colSums(buffel.comm) > 1] -> buffel.comm
  # transpose table
  buffel.comm.t <- t(buffel.comm)
  # co-occurrence analysis
  cooccur.buffel.sub <- cooccur(buffel.comm.t,
                          type = 'spp_site',
                          thresh = T,
                          spp_names = T)
  return(cooccur.buffel.sub)
}
  
# repeat function above 50 times and populate the results frame
for(i in 1:100){
  results.i <- random.introd(dataframe = buffel.int.data,
                             index.num = i)
  results.frame[i, 'total.pairs'] <- results.i$pairs
  results.frame[i, 'positive'] <- results.i$positive
  results.frame[i, 'positive.percent'] <- results.i$positive/results.i$pairs*100
  results.frame[i, 'negative'] <- results.i$negative
  results.frame[i, 'negative.percent'] <- results.i$negative/results.i$pairs*100
}

# mean species-pairs that occur more often than expected (positive associations)
pos <- as.array(results.frame$positive.percent)
pos.mean <- sum(pos)/nrow(results.frame)
pos.sd <- sd(pos)
  
# mean species-pairs that occur more often than expected (positive associations)
neg <- as.array(results.frame$negative.percent)
neg.mean <- sum(neg)/nrow(results.frame)
neg.sd <- sd(neg)

# write.csv(results.frame,
#           'SupplementaryTableS4.csv',
#           row.names = F)

# ---------------------------------------------------------------------#
## C: Fungal Trophic Mode ----
# ---------------------------------------------------------------------#

fungal_traits() -> fung.traits.db

### Read in taxonomic data
full.data <- read.csv('E13-133.data.output/fungal.taxonomy/TBAS/95sim_assignments_reportET66UU4R.csv')
full.data$otu <- tolower(full.data$otu)

### Isolate the C. ciliaris OTU only
# read in C. ciliaris dataset
ciliaris.data <- read.csv('E13-133.data.output/Range_OTU_Overlap.csv')

full.data[full.data$otu %in% ciliaris.data$OTU,] -> subset.data


### Compare out database to fung.traits.db ----
# database matches
fung.traits.db[fung.traits.db$Genus %in% subset.data$Genus.level.assignment,] -> match.db

# Study matches
subset.data[subset.data$Genus.level.assignment %in% fung.traits.db$Genus,] -> match.studies
colnames(match.studies)[colnames(match.studies) == "otu"] <- "OTU"

# Add information to ciliaris.data dataframe
left_join(ciliaris.data, match.studies, by = 'OTU') -> ciliaris.full

# Change those with NA or "" under Trophic.Mode to 'unclassified'
ciliaris.full$Trophic.Mode <- replace(ciliaris.full$Trophic.Mode,
                                      is.na(ciliaris.full$Trophic.Mode),
                                      "Unclassified")
ciliaris.full$Trophic.Mode <- replace(ciliaris.full$Trophic.Mode,
                                      ciliaris.full$Trophic.Mode == "",
                                      "Unclassified")

#### Create summary table for count of trophic level within native, introduced, and both ranges -----
# Unclassified Included
ciliaris.full %>%
  group_by(Range, Trophic.Mode) %>%
  summarise(count = n())  %>%
  group_by(Range) %>%
  mutate(percentage = count/sum(count) * 100) -> summary.range

# write.csv(summary.range,
#           'SupplementaryTableS3.csv',
#           row.names = F)

summary.range %>%
  filter(Trophic.Mode != 'Unclassified') %>%
  ggplot(aes(x = Range,
           y = count,
           fill = Trophic.Mode)) +
  geom_bar(stat = 'identity') +
  scale_fill_brewer(palette = 'BrBG') +
  xlab("") +
  ylab('Proposed trophic mode') +
  theme_bw() +
  theme(axis.text.y = element_text(size = 16,
                            color = 'black'),
        axis.text.x = element_text(size = 14,
                            color = 'black'),
        axis.title = element_text(size = 18, 
                                  color = 'black')) +
   labs(fill = "Trophic mode")

# ggsave('SupplementaryFigS2.tiff',
#        plot = last_plot(),
#        height = 5, 
#        width = 10, 
#        units = 'in')

#========================================================================================#
# Part 2: All data from introduced range ----------------
#========================================================================================#

otu.data <- read.csv('E13-133.data.output/SitexSpecies_95sim_Raref.csv')
clim.data <- as.data.frame(read_csv('E13-133.data.output/ClimateData.csv'))

#----------------------------------------------------------------------#
## D: Isolate data from introduced range  -----
#----------------------------------------------------------------------#

# Removed data from other plant parts than leaves (roots, dead leaves, litter)
# Removed data from plants that were washed or placed in ethanol

otu.data %>%
  filter(Location == 'US',
         Plant_Part == "Leaf Green",
         Wash == "Unwashed") -> us.data
clim.data[clim.data$Illumina %in% us.data$samples,] -> us.clim
us.data <- us.data[!colnames(us.data) %in% other.data$Otu]

# add metadata to us.data data.frame
for(i in us.data$samples){
  mat.i <- us.clim[us.clim$Illumina == i, 'BIO1']
  map.i <- us.clim[us.clim$Illumina == i, 'BIO12']
  bio10.i <- us.clim[us.clim$Illumina == i, 'BIO10']
  bio11.i <- us.clim[us.clim$Illumina == i, 'BIO11']
  bio16.i <- us.clim[us.clim$Illumina == i, 'BIO16']
  bio17.i <- us.clim[us.clim$Illumina == i, 'BIO17']
  lat.i <- us.clim[us.clim$Illumina == i, 'lat']
  long.i <- us.clim[us.clim$Illumina == i, 'long']
  host.i <- us.clim[us.clim$Illumina == i, 'Host.genus']
  us.data[us.data$samples == i, 'BIO1'] <- mat.i
  us.data[us.data$samples == i, 'BIO12'] <- map.i
  us.data[us.data$samples == i, 'BIO10'] <- bio10.i
  us.data[us.data$samples == i, 'BIO11'] <- bio11.i
  us.data[us.data$samples == i, 'BIO16'] <- bio16.i
  us.data[us.data$samples == i, 'BIO17'] <- bio17.i
  us.data[us.data$samples == i, 'lat'] <- lat.i
  us.data[us.data$samples == i, 'long'] <- long.i
  us.data[us.data$samples == i, 'Host.genus'] <- host.i
}

# reorder
us.data <- us.data[c(1:10,357:365,11:356)]

# Filter out columns with no occurrences
us.data[20:length(us.data)] -> comm.data
comm.data[colSums(comm.data) > 0] -> comm.data

# merge with metadata
cbind(us.data[1:19], comm.data) -> us.data

#----------------------------------------------------------------------#
## E: Dissimilarity comparison by host ----
#----------------------------------------------------------------------#

### Within and between group similarity ----
comm.data <- us.data[20: length(us.data)]
comm.data <- comm.data[colSums(comm.data) > 10]

## Overall
rownames(comm.data) <- us.data$samples # make sample names row names

# Jaccard dissimilarity
jacc.dist <- vegdist(comm.data, method = 'jaccard') 
jacc.dist.df <- as.data.frame(as.matrix(jacc.dist))# into data frame 
jacc.dist.df$samp1 <- rownames(jacc.dist.df) # make row names into column

# Make into long form
jacc.dist.df %>%
  gather(key = samp2, value = jaccard.dist, -samp1) -> overall.dist.long

# add host data
for(i in unique(overall.dist.long$samp1)){
  for(t in unique(overall.dist.long$samp2)){
    samp1.i <- us.data[us.data$samples == i, 'Host']
    samp2.i <- us.data[us.data$samples == t, 'Host']
    overall.dist.long[overall.dist.long$samp1 == i, 'Host1'] <- samp1.i
    overall.dist.long[overall.dist.long$samp2 == t, 'Host2'] <- samp2.i
  }
}

# Make third column for whether it is the same range or different
overall.dist.long[overall.dist.long$Host1 == overall.dist.long$Host2, 'Group.host'] <- 'Same'
overall.dist.long[overall.dist.long$Host1 != overall.dist.long$Host2, 'Group.host'] <- 'Different'

### Assess differences in similarity when grouped by native status -----
status <- data.frame(host = c("Cenchrus_ciliaris", "Bouteloua_curtipendula", "Hilaria_mutica", 
           "Pappophorum_bicolor", "Bothriochloa_ischaemum",
           "Heteropogon_contortus", "Sporobolus_cryptandrus", "Setaria",
           "Aristida", "Chloris", "Sporobolus", "Panicum_antidotale"),
           status = c("Nonnative", "Native", "Native",
                      "Native", "Nonnative",
                      "Nonnative", "Native", "Native",
                      "Native","Native","Native","Nonnative"))

for(i in 1:nrow(overall.dist.long)){
  host1.i <- overall.dist.long[i, 'Host1']
  status1.i <- status[status$host == host1.i, 'status']
  overall.dist.long[i, 'Status1'] <- status1.i
  host2.i <- overall.dist.long[i, 'Host2']
  status2.i <- status[status$host == host2.i, 'status']
  overall.dist.long[i, 'Status2'] <- status2.i
}

# Add Group data for native status
overall.dist.long[overall.dist.long$Status1 == overall.dist.long$Status2,
                  'Group.status'] <- 'Same'
overall.dist.long[overall.dist.long$Status1 != overall.dist.long$Status2,
                  'Group.status'] <- 'Different'

# Summary data of community similarity across native species, across non-native
# species, and across all species
overall.dist.long %>%
  filter(jaccard.dist > 0) %>%
  group_by(Status1, Status2, Group.status) %>%
  summarise(mean = mean(jaccard.dist),
            sd = sd(jaccard.dist)) -> jaccard.summary.native

# write.csv(jaccard.summary.native,
#           'IntroducedRange_NativeNonnativeDissimilarity.csv',
#           row.names = F)


# Kruskal-wallis: within native and within non-native
overall.dist.long %>%
  filter(jaccard.dist > 0,
         Group.status == 'Same') -> overall.dist.long.wtn.nat.nonnat

grp.glm <- kruskal.test(jaccard.dist ~ Status1,
                          data = overall.dist.long.wtn.nat.nonnat)

# plot
overall.dist.long.wtn.nat.nonnat %>%
  filter(jaccard.dist > 0) %>%
  ggplot(aes(x = Status1,
             y = jaccard.dist,
             fill = Status1)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  geom_point(size = 0.5) +
  theme_classic() +
  xlab('') +
  ylab('Jaccard dissimilarity') +
  facet_grid(. ~ Status2) +
  scale_fill_brewer(palette = "BuGn") +
  theme(text = element_text(size=20, color = 'black'),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank())

# ggsave('Figure6_WithinGroup.tiff',
#        plot = last_plot(), device = 'tiff', scale = 2,
#        width = 14, height = 10, unit = 'cm')

# Kruskal-wallis: within same group and between groups
# factor explanatory variables
overall.dist.long$Status1 <- as.factor(overall.dist.long$Status1)
overall.dist.long$Group.status <- as.factor(overall.dist.long$Group.status)

overall.dist.long %>%
  filter(jaccard.dist > 0) -> overall.dist.long.wtn.btwn

grp.glm <- kruskal.test(jaccard.dist ~ Group.status,
                          data = overall.dist.long.wtn.btwn)

# plot
overall.dist.long.wtn.btwn %>%
  filter(jaccard.dist > 0,
         Group.status == 'Different') %>%
  ggplot(aes(x = Group.status,
             y = jaccard.dist)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  geom_point(size = 0.5) +
  theme_classic() +
  xlab('') +
  ylab('Jaccard dissimilarity') +
  # facet_grid(. ~ Status2) +
  scale_fill_brewer(palette = "BuGn") +
  theme(text = element_text(size=20, color = 'black'),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank())

# ggsave('Figure6_BetweenGroup.tiff',
#        plot = last_plot(), device = 'tiff', scale = 2,
#        width = 4, height = 10, unit = 'cm')

#----------------------------------------------------------------------#
## F: Taxonomic groups shared by all grass species in introduced range ----
#----------------------------------------------------------------------#

# Pull out C. ciliaris OTU to compare
us.data %>%
  select(Host, starts_with('Otu')) %>%
  gather(key = Otu, value = read.ab, -Host) %>%
  group_by(Host, Otu) %>%
  summarise(read.ab = sum(read.ab)) %>%
  filter(read.ab > 0, 
         Host == "Cenchrus_ciliaris") %>%
  select(Otu) -> cenchrus.us.otu

# Pull out OTU associated with co-occurring grasses in the introduced range
us.data %>%
  select(Host, starts_with('Otu')) %>%
  gather(key = Otu, value = read.ab, -Host) %>%
  group_by(Host, Otu) %>%
  summarise(read.ab = sum(read.ab)) %>%
  filter(read.ab > 0, 
         Host != "Cenchrus_ciliaris") -> grass.us.otu

unique(grass.us.otu$Otu) -> grass.us.unique.otu

# pull out otu shared with Cenchrus and other grass epcies in the introduced range
cenchrus.us.otu[cenchrus.us.otu$Otu %in% grass.us.unique.otu,] -> shared.otu

# Read in taxonomic data
full.data <- read.csv('E13-133.data.output/fungal.taxonomy/TBAS/95sim_assignments_reportET66UU4R.csv')

# Create proportion data
# change OTU column to lower case in both full.data and summary.overlap files
full.data$otu <- tolower(full.data$otu)
shared.otu$Otu <- tolower(shared.otu$Otu)

for(i in shared.otu$Otu){
    class.i <- unique(full.data[full.data$otu == i, 'Most.common.class.level.assignment'])
    order.i <- unique(full.data[full.data$otu == i, 'Most.common.order.level.assignment'])
    family.i <- unique(full.data[full.data$otu == i, 'Most.common.family.level.assignment'])
    genus.i <- unique(full.data[full.data$otu == i, 'Most.common.Genus.level.assignment'])
    shared.otu[shared.otu$Otu == i, 'Class'] <- class.i
    shared.otu[shared.otu$Otu == i, 'Order'] <- order.i
    shared.otu[shared.otu$Otu == i, 'Family'] <- family.i
    shared.otu[shared.otu$Otu == i, 'Genus'] <- genus.i[1]
}

# Change order of classes
shared.otu %>%
  filter(Class %in% c("Dothideomycetes","Eurotiomycetes",
                      "Orbiliomycetes","Leotiomycetes",
                      "Pezizomycetes","Saccharomycetes",
                      "Sordariomycetes","Xylonomycetes")) -> shared.otu

shared.otu$Class <- factor(shared.otu$Class,
                           levels = c("Dothideomycetes","Eurotiomycetes",
                                      "Orbiliomycetes","Leotiomycetes",
                                      "Pezizomycetes","Saccharomycetes",
                                      "Sordariomycetes","Xylonomycetes"))

# Create a matrix with counts
shared.otu %>%
  group_by(Class) %>%
  summarise(n = n()) -> shared.otu.summary

shared.otu.summary$Group <- 'Shared_IntroducedRange'

### GGplot Percent stacked: Shared in introduced range across grasses ----
class.plot <- shared.otu.summary %>%
  group_by(Group) %>%
  mutate(RelativeAbundance = n / sum(n)) %>%
   ggplot(aes(x = Group,
             y = RelativeAbundance,
             fill = Class)) +
  geom_bar(stat = 'identity') +
  scale_fill_brewer(palette = 'RdBu') +
  xlab('') +
  ylab('Relative abundance of OTUs') +
  theme_bw() +
  theme(axis.text.y = element_text(size = 16,
                            color = 'black'),
        axis.title = element_text(size = 18, 
                                  color = 'black'),
        axis.text.x = element_blank(),
        legend.position = 'none')

# ggsave('Figure3B.tiff',
#        plot = class.plot,
#        width = 7.5,
#        height = 20,
#        units = 'cm')

### Chi-square: comparing differences in taxonomic composition of Cenchrus to native and non-native grasses
### 
### [p]
us.data %>%
  select(Host, Range, starts_with('Otu')) %>%
  gather(key = Otu, value = read.ab, -Host, -Range) %>%
  group_by(Host, Range, Otu) %>%
  summarise(read.ab = sum(read.ab)) %>%
  filter(read.ab > 0) -> us.summary

us.summary$Otu <- tolower(us.summary$Otu)

# add taxonomic data
for(i in us.summary$Otu){
    class.i <- unique(full.data[full.data$otu == i, 'Most.common.class.level.assignment'])
    order.i <- unique(full.data[full.data$otu == i, 'Most.common.order.level.assignment'])
    family.i <- unique(full.data[full.data$otu == i, 'Most.common.family.level.assignment'])
    genus.i <- unique(full.data[full.data$otu == i, 'Most.common.Genus.level.assignment'])
    us.summary[us.summary$Otu == i, 'Class'] <- class.i
    us.summary[us.summary$Otu == i, 'Order'] <- order.i
    us.summary[us.summary$Otu == i, 'Family'] <- family.i
    us.summary[us.summary$Otu == i, 'Genus'] <- genus.i[1]
}

# summary of OTU by class and host
us.summary %>%
  group_by(Host, Range, Class) %>%
  summarise(count = n()) -> us.summary.chi

# create a new column Cenchrus and non-Cenchrus
for(i in 1:nrow(us.summary.chi)){
  if(us.summary.chi[i, 'Host'] != 'Cenchrus_ciliaris'){
    us.summary.chi[i, 'Host.o'] <- 'Non-cenchrus' 
  } else(us.summary.chi[i, 'Host.o'] <- 'Cenchrus' )
}

# Create native, invasive, and cenchrus column
for(i in 1:nrow(us.summary.chi)){
  if(us.summary.chi[i, 'Range'] == 'Invasive' & us.summary.chi[i, 'Host'] != 'Cenchrus_ciliaris'){
    us.summary.chi[i, 'Range.o'] <- 'Invasive'}
  if(us.summary.chi[i, 'Range'] == 'Invasive' & us.summary.chi[i, 'Host'] == 'Cenchrus_ciliaris'){
    us.summary.chi[i, 'Range.o'] <- 'Cenchrus'} 
  if(us.summary.chi[i, 'Range'] != 'Invasive'){
    us.summary.chi[i, 'Range.o'] <- 'Native'
  }
}

# Group by new column
us.summary.chi %>%
  group_by(Range.o, Class) %>%
  summarise(count.total = sum(count)) %>%
  filter(Class %in% c("Dothideomycetes","Eurotiomycetes",
                      "Orbiliomycetes","Leotiomycetes",
                      "Pezizomycetes","Saccharomycetes",
                      "Sordariomycetes","Xylonomycetes")) %>%
  spread(key = Class, value = count.total, fill = 0) -> us.summary.chi.broad

us.summary.chi.broad <- as.data.frame(us.summary.chi.broad)
row.names(us.summary.chi.broad) <- us.summary.chi.broad$Range.o
us.summary.chi.broad[-1] -> class.chi

chisq.test(class.chi)

#----------------------------------------------------------------------#
## G: Comparison of OTU shared between ranges with OTU found in introduced range ----
#----------------------------------------------------------------------#

# read in data frame with OTU shared between ranges and filter out OTU unique
# to native or introduced range, only keep shared OTU
shared.data <- read.csv('E13-133.data.output/Range_OTU_Overlap.csv')
shared.data <- filter(shared.data, Range == 'Shared')

# remove C. ciliaris data and group data by host species
us.data %>%
  filter(Host != 'Cenchrus_ciliaris') %>%
  select(Host, Range, starts_with('Otu')) %>%
  group_by(Host, Range) %>%
  summarise(across(starts_with('Otu'), ~sum(., na.rm = TRUE))) -> compare.data

# remove columns with no reads
compare.data[3:length(compare.data)] -> clean.data
clean.data[colSums(clean.data) > 0] -> clean.data
cbind(compare.data[c('Host', 'Range')], clean.data) -> compare.data

# create data frame with absence data removed so it can be compared to shared.data
compare.data %>%
  gather(key = OTU, value = n, -Host, -Range) %>%
  filter(n > 0) -> compare.data

# make all OTU names lower case
compare.data$OTU <- tolower(compare.data$OTU)

### which OTU from the US species is shared between Kenya and US ------
compare.data[compare.data$OTU %in% shared.data$OTU,] -> compare.shared

# write.csv(compare.shared,
#           'IntroducedRange_NonCiliaris_Shared.csv',
#           row.names = F)

#### Assess how many of these were found on ----
compare.shared %>%
  group_by(OTU) %>%
  summarise(number.hosts = n(),
            total.reads = sum(n)) -> summary

compare.shared %>%
  group_by(Range, OTU) %>%
  summarise(number.hosts = n(),
            total.reads = sum(n)) %>%
  ggplot(aes(x = OTU,
           y = number.hosts,
           fill = Range)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = c('#636363','#3182bd')) +
  xlab("") +
  ylab('Number of Non-ciliaris hosts') +
  theme_bw() +
  theme(axis.text.y = element_text(size = 16,
                            color = 'black'),
        axis.text.x = element_text(size = 10,
                            color = 'black',
                            angle = 45,
                            hjust = 1.0),
        axis.title = element_text(size = 18, 
                                  color = 'black')) +
   labs(fill = "Host native status")

ggsave('SupplementaryFigS3.tiff',
       plot = last_plot(), 
       width = 30, 
       height = 12, 
       units = 'cm')

### which OTU from the shared species is in the US ----
shared.data[shared.data$OTU %in% compare.data$OTU,] -> shared.shared
shared.shared %>%
  group_by(Class) %>%
  summarise(count = n())

### which OTU from the shared species is only on C. ciliaris -----
shared.data[!shared.data$OTU %in% compare.data$OTU,] -> shared.Notshared
shared.Notshared %>%
  group_by(Class) %>%
  summarise(count = n())

# write.csv(shared.Notshared,
#           'IntroducedRange_UniqueToCiliaris.csv',
#           row.names = F)

# ----------------------------------------------------------------------#
## Summary stats ----
# ----------------------------------------------------------------------#

## OTUs per species, shared and unique: all otu ----
us.data %>%
  dplyr::select(-samples, -Latitude, -Longitude, -Year, -Site, -Wash,
                -Plant_Part, -Location, -Range, -BIO1, -BIO12, -BIO10,
                -BIO11, -BIO16, -BIO17, -lat, -long, -Host.genus) %>%
  gather(key = 'OTU', value = 'read.number', -Host) %>%
  filter(read.number > 0) -> us.filt

status <- data.frame(host = c("Cenchrus_ciliaris", "Bouteloua_curtipendula", "Hilaria_mutica", 
           "Pappophorum_bicolor", "Bothriochloa_ischaemum",
           "Heteropogon_contortus", "Sporobolus_cryptandrus", "Setaria",
           "Aristida", "Chloris", "Sporobolus", "Panicum_antidotale"),
           status = c("Introduced", "Native", "Native",
                      "Native", "Nonnative",
                      "Nonnative", "Native", "Native",
                      "Native","Native","Native","Nonnative"))

# add native status to us.filt
for(h in status$host){
  us.filt[us.filt$Host == h, 'status'] <- status[status$host == h, 'status']
}

# isolate OTU unique to each group (Cenchrus cilaiaris)
us.filt %>%
  filter(Host == 'Cenchrus_ciliaris') %>%
  distinct(OTU)-> us.buf
us.filt %>%
  filter(status == 'Native') %>%
  distinct(OTU) -> us.nat
us.filt %>%
  filter(status == 'Nonnative') %>%
  distinct(OTU) -> us.nonnat

# Buffel compared to native and nonnative
# Shared with native 
us.buf[us.buf$OTU %in% us.nat$OTU,] -> shared.buf.nat

# Shared with non-native 
us.buf[us.buf$OTU %in% us.nonnat$OTU,] -> shared.buf.nonnat

#Unique
us.buf[!us.buf$OTU %in% us.nat$OTU & !us.buf$OTU %in% us.nonnat,] -> uniquetobuff.natnonnat

# Unique to native
us.nat[!us.nat$OTU %in% us.buf$OTU & !us.nat$OTU %in% us.nonnat$OTU,] -> uniquetonat.bufnonnat

# Unique to non-native
us.nonnat[!us.nonnat$OTU %in% us.buf$OTU & !us.nonnat$OTU %in% us.nat$OTU,] -> uniquetononnat.bufnat

# Compare the OTU from C. cilairis shared between ranges to those found in the introduced range ----
range.shared <- read.csv('E13-133.data.output/Range_OTU_Overlap.csv')

# isolate shared otu
range.shared <- filter(range.shared, Range == 'Shared')

# isolate all Otus not within Cenchrus ciliaris
filter(us.data, Host.genus != 'Cenchrus') -> us.noncenchrus.data

# make sample name row name
rownames(us.noncenchrus.data) <- us.noncenchrus.data$samples

# remove meta data
us.otus <- us.noncenchrus.data[23:length(us.noncenchrus.data)]
us.otus.meta <- us.noncenchrus.data[1:22]

# remove otus with no data
us.otus[colSums(us.otus) > 0] -> us.otus

# isolate otu names
tolower(names(us.otus)) -> us.otu.names

## Otus found in other native and non-native in introduced range ----
range.shared[range.shared$OTU %in% us.otu.names,] -> shared.introduced

## Otus nont found in other native and non-native in introduced range ----
range.shared[!range.shared$OTU %in% us.otu.names,] -> unique.introduced

# OTUs at the order level -----
# create dataframe with order data
us.data %>%
  select(-Wash, -Plant_Part, -Location, -Latitude, -Longitude, -Year, -samples) %>%
  gather(key = Query.sequence, value = reads, -Host, -Site, -Range) %>%
  group_by(Query.sequence, Host, Range) %>%
  summarise(total = sum(reads)) %>%
  filter(total > 0) -> order.data

left_join(order.data, full.data) -> order.data.b

# Create count table of orders by host species
order.data.b %>%
  group_by(Host, Most.common.order.level.assignment) %>%
  summarise(count = n()) -> order.summary

order.summary %>%
  group_by(Host) %>%
  summarise(total_count = sum(count)) -> host_total

order.summary %>%
  left_join(host_total, by = "Host") %>%
  mutate(proportion = count/total_count) %>%
  filter(Most.common.order.level.assignment %in%
           c("Pleosporales", "Capnodiales", "Eurotiales", "Hypocreales")) -> order.summary.b

# Plot
ggplot(order.summary.b,
       aes(x = Host,
           y = proportion,
           fill = Most.common.order.level.assignment)) +
  geom_bar(stat = "identity") +
  xlab("") +
  ylab("Proportion") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black",
                                   angle = 45, hjust = 1))

# average of order proportion across all host species
order.summary.b %>%
  group_by(Most.common.order.level.assignment) %>%
  summarise(mean.proportion = mean(proportion),
            sd.proportion = sd(proportion))

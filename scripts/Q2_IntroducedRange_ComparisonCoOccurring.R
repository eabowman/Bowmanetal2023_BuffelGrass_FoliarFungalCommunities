## Script created by Dr. Liz Bowman (eabowman@utexas.edu), May 17, 2021 
## Focus is to answer question 2: what drives foliar fungal community
## composition within the introduced range compared to co-occurring native
## and non-native plants
## 
## Data: 95% Sequence Similarity Rarefied data (SitexSpecies_95sim_Raref.csv)
## A: Isolate data from hosts in the introduced range
## B: Calculate phylogenetic distance of host plants in introduced range and
## geographic distance between sites
## C: PERMANOVA with precipitation, temperature, geographic distance, and
## phylogenetic distance

# Be sure to run the LoadLibraries.R file before running code below.

otu.data <- read.csv('data/SitexSpecies_95sim_Raref.csv')
clim.data <- as.data.frame(read_csv('data/ClimateData.csv'))
other.data <- read_csv('data/eukaryote_zotu95-Phylum.csv')

#----------------------------------------------------------------------#
# A: Isolate data from introduced range  -----
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
# B: Calculate phylogenetic and geographic distance between sites>> -----
#----------------------------------------------------------------------#

## Phylogenetic distance-----

## Create distance matrix
matK.tree <- read.tree('data/trnKmatK_branchlengths_bootstrapValues.nwk')
dist.matrix.2d <-cophenetic.phylo(matK.tree)
# Principal coordinates analysis (PCoA) to transform the pairwise patristic 
# distance matrix to phylogenetic eigenvectors 
host.pcoa <- cmdscale(dist.matrix.2d)
# change to genus instead of species
rownames(host.pcoa) <- c('Cenchrus', 'Setaria', 'Panicum', 'Bothriochloa', 
                         'Heteropogon', 'Aristida', 'Sporobolus', 'Chloris',
                         'Pappophorum', 'Bouteloua', 'Hilaria')

# add values to  file so it matches size
for(i in rownames(host.pcoa)){
  v1.i <- host.pcoa[rownames(host.pcoa) == i, 1]
  v2.i <- host.pcoa[rownames(host.pcoa) == i, 2]
  us.data[us.data$Host.genus == i, 'host.pcoa.1'] <- v1.i
  us.data[us.data$Host.genus == i, 'host.pcoa.2'] <- v2.i
}

## Compute distance matrix for environmental data-----
# Historical climate: MAP, MAT, MPWQ, MTWQ
us.data %>%
  filter(Location == 'US') %>%
  rename(MAP = BIO12, MAT = BIO1, MPWQ = BIO16, MTWQ = BIO10) %>%
  dplyr::select(samples, MAP, MAT, MPWQ, MTWQ) -> clim.meta.hist
rownames(clim.meta.hist) <- clim.meta.hist$samples
clim.hist.dist <- vegdist(clim.meta.hist[2:length(clim.meta.hist)],
                          method = 'euclidean')

# Recent climate: Tmax.18, Prec.18
us.clim %>%
  filter(Location == 'US') %>%
  dplyr::select(Illumina, Tmax.18, Prec.18) -> clim.meta.rec
rownames(clim.meta.rec) <- clim.meta.rec$Illumina
clim.rec.dist <- vegdist(clim.meta.rec[2:length(clim.meta.rec)],
                         method = 'euclidean')

## Compute PCNM eigenvectors for location data----
# Geographic distance
geo.dist <- vegdist(us.data[c('lat','long')], method = 'euclidean')
geo.pcnm <- pcnm(geo.dist)

# add to us.clim data frame
us.data$geo.pcnm.1 <- geo.pcnm$vectors[,1]

# reorder
names(us.data)
us.data <- us.data[c(1:19,245:247,20:244)]

#----------------------------------------------------------------------#
# C: PERMANOVA with precipitation, temperature, geographic distance, and phylogenetic distance >> ----
#----------------------------------------------------------------------#

## Assess correlation of native status with phylogenetic distance ----
## Based on the phylogeny, these two factors seem to be confounded.

kruskal.test(host.pcoa.1 ~ Range, data = us.data)

## PERMANOVA ----
# Isolate data
comm.data <- us.data[23: length(us.data)]
comm.data <- comm.data[colSums(comm.data) > 10]

# Hellinger transformation
comm.hell <- decostand(comm.data, method = 'total')

# Jaccard dissimilarity
jacc.dist <- vegdist(comm.hell, method = 'jaccard') 

# Forward selection of explanatory variables
# Create data frame with all possible explanatory variables
us.data %>%
  select(BIO1, BIO12, BIO16, BIO17, BIO10, BIO11, host.pcoa.1, geo.pcnm.1) -> expl.data
forward.sel(comm.data, expl.data)

# PERMANOVA: Table 4
jacc.adonis <- adonis2(jacc.dist ~ host.pcoa.1 * BIO12 * geo.pcnm.1,
                      data = us.data,
                      strata = us.data$Site)
jacc.adonis

# write results
jacc.adonis.results <- as.data.frame(jacc.adonis)
# write.csv(jacc.adonis.results,
#           'Table2.csv',
#           row.names = T)

## NMDS plot
jacc.mds <- metaMDS(jacc.dist, dist = 'bray',
                    try = 1000, trymax = 1000,
                    tidy = T, shrink = T, k = 2)
jacc.stress <- jacc.mds$stress

# format data for plot
us.data$Host.genus <- factor(us.data$Host.genus,
                       levels = c("Cenchrus", "Setaria", "Panicum", "Bothriochloa",
                                  "Heteropogon", "Aristida", "Sporobolus", "Chloris",
                                  "Pappophorum", "Bouteloua", "Hilaria"))
data.scores <- data.frame(NMDS1 = jacc.mds$points[,1],
                              NMDS2 = jacc.mds$points[,2],
                              site = us.data$Site,
                              Range = us.data$Range,
                              Host = us.data$Host.genus, 
                              MAP = us.data$BIO12)


# Figure 4A: By MAP
jacc.plot <- ggplot() +
  geom_point(data = data.scores, aes(x = NMDS1,
                                     y = NMDS2,
                                     color = MAP),
             size = 3) +
  scale_colour_gradient(low = "black", high = "lightblue") +
  coord_equal() +
  theme_classic() +
  theme(axis.text = element_text(size = 14,
                                 colour = 'Black'),
        axis.title = element_text(size = 16))

# ggsave('Fig4A.tiff',
#        plot = jacc.plot, width = 20, height = 20, units = 'cm')


# Figure 4B: By host
jacc.plot <- ggplot() +
  geom_point(data = data.scores, aes(x = NMDS1,
                                     y = NMDS2,
                                     color = Host),
             size = 3) +
  scale_color_brewer(palette = 'RdYlGn') +
  coord_equal() +
  theme_classic() +
  theme(axis.text = element_text(size = 14,
                                 colour = 'Black'),
        axis.title = element_text(size = 16))

# ggsave('Fig4B.tiff',
#        plot = jacc.plot, width = 20, height = 20, units = 'cm')
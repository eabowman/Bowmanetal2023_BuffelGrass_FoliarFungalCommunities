## Script created by Dr. Liz Bowman (eabowman@utexas.edu), May 14, 2021 
## Focus is to answer question 1: how an invasive plant’s foliar fungal 
## community overlaps between its native and introduced ranges. This includes:
## a) how similar OTU richness and community composition are between ranges and
## b) do similar drives shape community composition between the two ranges
## 
## Data: 95% Sequence Similarity Rarefied data (SitexSpecies_95sim_Raref.csv)
## A: Isolate data from C. ciliaris hosts in TX, AZ, and Kenya
## B: Species richness
## C: Community turnover and major taxonomic groups
## D: Environmental factors driving community composition
## E: Dissimilarity plots

# Be sure to run the LoadLibraries.R file before running code below.

# read in data
otu.data <- read.csv('data/SitexSpecies_95sim_Raref.csv')
clim.data <- read_csv('data/ClimateData.csv')
other.data <- read_csv('data/eukaryote_zotu95-Phylum.csv')

#========================================================================================#
# A: Isolate data from C. ciliaris hosts in TX, AZ, and Kenya----
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
# B: Species richness----
#========================================================================================#

## Species richness by range ====
# isolate community data
comm.data <- buffel.data[22:length(buffel.data)]
comm.data <- comm.data[colSums(comm.data) > 10]

# Make matrix with species richness
buffel.sr <- buffel.data[1:21]
buffel.sr$spec.richness <- specnumber(comm.data)

# assess normality
hist(buffel.sr$spec.richness)
shapiro.test(buffel.sr$spec.richness) # normal

### Anova with site nested as a random variable ----
# Mixed effect model
lme.sp <- lme4::lmer(spec.richness ~ Range + (1|Site), data = buffel.sr)
summary(lme.sp)

qqnorm(resid(lme.sp))
qqline(resid(lme.sp))

anova(lme.sp)

## Average species richness within each range ----

summary(buffel.sr$spec.richness)
boxplot(buffel.sr$spec.richness)

buffel.sr %>%
  group_by(Range) %>%
  summarise(mean = mean(spec.richness),
            median = median(spec.richness),
            sd = sd(spec.richness),
            min = min(spec.richness),
            max = max(spec.richness))

#========================================================================================#
# C: Community turnover and major taxonomic groups----
#========================================================================================#

## Community turnover ----
buffel.data %>%
  dplyr::select(Range, starts_with('Otu')) %>%
  gather(key = 'OTU', value = 'read.abund', -Range) %>%
  filter(read.abund > 0) %>%
  distinct(Range, OTU, .keep_all = T) -> buffel.sum

# Shared and Unique to native and Invasive ranges
buffel.sum %>%
  filter(Range == 'Native') %>%
  group_by(OTU) -> buffel.nat

buffel.sum %>%
  filter(Range == 'Invasive') %>%
  group_by(OTU) -> buffel.inv

buffel.inv[buffel.inv$OTU %in% buffel.nat$OTU, ] -> shared.buff
buffel.inv[!buffel.inv$OTU %in% buffel.nat$OTU, ] -> unique.inv.buff
buffel.nat[!buffel.nat$OTU %in% buffel.inv$OTU, ] -> unique.nat.buff

# Create Otu table
shared.buff$Range <- 'Shared'
unique.inv.buff$Range <- 'Introduced'
unique.nat.buff$Range <- 'Native'
rbind(shared.buff, unique.inv.buff, unique.nat.buff) -> summary.overlap

# read in taxonomic data
full.data <- read.csv('data/95sim_assignments_reportET66UU4R.csv')
# change OTU column to lower case in both full.data and summary.overlap files
full.data$otu <- tolower(full.data$otu)
summary.overlap$OTU <- tolower(summary.overlap$OTU)

for(i in summary.overlap$OTU){
  if(i %in% full.data$otu){
    class.i <- unique(full.data[full.data$otu == i, 'Most.common.class.level.assignment'])
    order.i <- unique(full.data[full.data$otu == i, 'Most.common.order.level.assignment'])
    family.i <- unique(full.data[full.data$otu == i, 'Most.common.family.level.assignment'])
    genus.i <- unique(full.data[full.data$otu == i, 'Most.common.Genus.level.assignment'])
    summary.overlap[summary.overlap$OTU == i, 'Class'] <- class.i
    summary.overlap[summary.overlap$OTU == i, 'Order'] <- order.i
    summary.overlap[summary.overlap$OTU == i, 'Family'] <- family.i
    summary.overlap[summary.overlap$OTU == i, 'Genus'] <- genus.i[1]
  }
}

# write.csv(summary.overlap, 'Range_OTU_Overlap.csv', row.names = F)

# Total per range (Introduced, native, shared)
summary.overlap %>%
  group_by(Range) %>%
  count() -> range.count

# Class data
summary.overlap %>%
  group_by(Range,Class) %>%
  count()  -> class.count

for(r in 1:nrow(class.count)){
  if(class.count[r, 'Range'] == 'Introduced'){
    class.count[r, 'rel.ab'] <- class.count[r, 'n']/108
  }
  if(class.count[r, 'Range'] == 'Native'){
    class.count[r, 'rel.ab'] <- class.count[r, 'n']/52
  }
  if(class.count[r, 'Range'] == 'Shared'){
    class.count[r, 'rel.ab'] <- class.count[r, 'n']/76
  }
}

# write.csv(class.count,
#           'Range_OTU_Overlap_ClassSummary.csv',
#           row.names = F)

## Taxonomic data: class level ----

class.count %>%
  filter(n > 1) %>%
  filter(Class %in% c('Dothideomycetes', 'Eurotiomycetes',
                      'Orbiliomycetes', 'Leotiomycetes',
                      'Pezizomycetes', 'Saccharomycetes',
                      'Sordariomycetes', 'Xylonomycetes')) %>%
  # Check Lecanoromycetes, Neolectomycetes, Schizosaccharomycetes
   group_by(Range) %>%
   mutate(RelativeAbundance = n / sum(n)) %>%
   ggplot(aes(x = Range,
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
        legend.position = 'right')

# ggsave('Fig3A.tiff',
#        plot = last_plot(),
#        width = 20,
#        height = 20,
#        units = 'cm')

### Chi-square ----
# Remove shared data
class.count %>%
  # filter(Class %in% c('Dothideomycetes', 'Eurotiomycetes',
  #                     'Orbiliomycetes', 'Leotiomycetes',
  #                     'Pezizomycetes', 'Saccharomycetes',
  #                     'Sordariomycetes', 'Xylonomycetes')) %>% # same classes as graph
  filter(Range != 'Shared') %>%
  select(-rel.ab) %>%
  spread(key = Class, value = n, fill = 0) -> class.chi
class.chi <- as.data.frame(class.chi)
row.names(class.chi) <- class.chi$group
class.chi[-1] -> class.chi

chisq.test(class.chi)

#========================================================================================#
# D. Environmental factors driving community composition ----
#========================================================================================#

## PERMANOVA with precipitation, temperature, geographic location, year----

# Remove outliers KMBGFF1
buffel.data.out <- buffel.data[!buffel.data$samples == 'KMBGFF1',]

# isolate community data and create jaccard distance matrix
comm.data <- buffel.data.out[22:length(buffel.data.out)]
comm.data <- comm.data[colSums(comm.data) > 10]

env.permanova <- adonis2(comm.data ~ Range * Geo.pcnm.2 * MAP, 
                        method = 'jaccard',
                        strata = buffel.data.out$Site,
                        data = buffel.data.out)
env.permanova

# write results
jacc.adonis.results <- as.data.frame(env.permanova)
# write.csv(jacc.adonis.results,
#           'results/Table1.csv',
#           row.names = T)

## NMDS plot
jacc.dist <- vegdist(comm.data, method = 'jaccard', binary = T)
jacc.mds <- metaMDS(jacc.dist, dist = 'bray',
                    try = 1000, trymax = 1000)
jacc.stress <- jacc.mds$stress

# format data for plot
data.scores <- data.frame(NMDS1 = jacc.mds$points[,1],
                          NMDS2 = jacc.mds$points[,2],
                          state = buffel.data.out$Location.state,
                          site = buffel.data.out$Site,
                          Range = buffel.data.out$Range,
                          MAP = buffel.data.out$MAP,
                          MAT = buffel.data.out$MAT)

jacc.plot <- ggplot() + 
  geom_point(data = data.scores, aes(x = NMDS1,
                                     y = NMDS2,
                                     shape = Range,
                                     color = state),
             size = 4) +
  scale_color_manual(values = c('#9ecae1','#bdbdbd','#636363','#3182bd')) +
  coord_equal() +
  theme_classic() +
  theme(axis.text = element_text(size = 14,
                                 colour = 'Black'),
        axis.title = element_text(size = 16))

# ggsave('figures/Fig2A.tiff',
#        plot = jacc.plot, width = 20, height = 20, units = 'cm')

#========================================================================================#
# E: Dissimilarity plots----
#========================================================================================#

# isolate community data
comm.data <- buffel.data[22:length(buffel.data)]

# Remove OTU with less than 10 occurrences
comm.data <- comm.data[colSums(comm.data) > 10]

# create distance matrix
jacc.dist <- vegdist(comm.data, method = "jaccard", binary = T)

# assess homogeneity of variance: NOT NORMAL
jacc.betadisper <- betadisper(jacc.dist, group = buffel.data$Location)
anova(jacc.betadisper)

## ANOSIM: Significant, explains ~47.1% of variation -----
range.anosim <- anosim(jacc.dist, grouping = buffel.data$Range,
                          distance = 'jaccard')
range.anosim

## Within and between group similarity ----
## Overall
# isolate community data
comm.data <- buffel.data[22:length(buffel.data)]
comm.data <- comm.data[colSums(comm.data) > 10]
rownames(comm.data) <- buffel.data$samples # make sample names row names

# Jaccard dissimilarity
jacc.dist <- vegdist(comm.data, method = 'jaccard') 
jacc.dist.df <- as.data.frame(as.matrix(jacc.dist))# into data frame 
jacc.dist.df$samp1 <- rownames(jacc.dist.df) # make row names into column

# Make into long form
jacc.dist.df %>%
  gather(key = samp2, value = jaccard.dist, -samp1) %>%
  filter(jaccard.dist > 0) -> overall.dist.long

# add range data
for(i in unique(overall.dist.long$samp1)){
  for(t in unique(overall.dist.long$samp2)){
    samp1.i <- clim.data[clim.data$Illumina == i, 'Range']
    samp2.i <- clim.data[clim.data$Illumina == t, 'Range']
    overall.dist.long[overall.dist.long$samp1 == i, 'Range1'] <- samp1.i
    overall.dist.long[overall.dist.long$samp2 == t, 'Range2'] <- samp2.i
    site1.i <- clim.data[clim.data$Illumina == i, 'Site']
    site2.i <- clim.data[clim.data$Illumina == i, 'Site', ]
    overall.dist.long[overall.dist.long$samp1 == i, 'Site1'] <- site1.i
    overall.dist.long[overall.dist.long$samp2 == i, 'Site2'] <- site2.i
    mpwq1.i <- clim.data[clim.data$Illumina == i, 'BIO16']
    mpwq2.i <- clim.data[clim.data$Illumina == i, 'BIO16']
    overall.dist.long[overall.dist.long$samp1 == i, 'MPWQ.1'] <- mpwq1.i
    overall.dist.long[overall.dist.long$samp2 == i, 'MPWQ.2'] <- mpwq2.i
    map1.i <- clim.data[clim.data$Illumina == i, 'BIO12']
    map2.i <- clim.data[clim.data$Illumina == i, 'BIO12']
    overall.dist.long[overall.dist.long$samp1 == i, 'MAP.1'] <- map1.i
    overall.dist.long[overall.dist.long$samp2 == i, 'MAP.2'] <- map2.i
  }
}

# Make third column for whether it is the same range or different
overall.dist.long[overall.dist.long$Range1 == overall.dist.long$Range2, 'Group'] <- 'Same'
overall.dist.long[overall.dist.long$Range1 != overall.dist.long$Range2, 'Group'] <- 'Different'

# Make third column for whether it is the same site or different
overall.dist.long[overall.dist.long$Site1 == overall.dist.long$Site2, 'Group.site'] <- 'Same'
overall.dist.long[overall.dist.long$Site1 != overall.dist.long$Site2, 'Group.site'] <- 'Different'
overall.dist.long -> site.dist.long

overall.dist.long <- filter(overall.dist.long, jaccard.dist > 0.4)

# assess normality
histogram(overall.dist.long$jaccard.dist) # not normal, strong left skew

### Assess similarity within and between ranges ----

# Kruskal-Wallis test
grp.wil <- kruskal.test(jaccard.dist ~ Group,
                       data = overall.dist.long)
grp.wil

# Means of each site
overall.dist.long %>%
  group_by(Group, Range1) %>%
  summarise(mean = mean(jaccard.dist),
            sd = sd(jaccard.dist)) -> jaccard.summary

# write.csv(jaccard.summary,
#           'NativeIntrodRange_CenchrusCiliarisDissimilarity.csv',
#           row.names = F)

## GGplot
# Between the two ranges
btwn.overall.dist.long <- filter(overall.dist.long, Group == 'Different')
btwn.plot <- ggplot(btwn.overall.dist.long,
                     aes(x = Group,
                         y = jaccard.dist)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        text = element_text(size=20),
        legend.position = "none")

# Within the each range type (Kenya and US)
wtn.overall.dist.long <- filter(overall.dist.long, Group == 'Same')

# Wilcoxon test
wtn.overall.dist.long[wtn.overall.dist.long$Range1 == 'Invasive', 'Code'] <- 1
wtn.overall.dist.long[wtn.overall.dist.long$Range1 == 'Native', 'Code'] <- 0
kruskal.test(wtn.overall.dist.long$jaccard.dist, wtn.overall.dist.long$Code)

wtn.plot <- ggplot(wtn.overall.dist.long,
                     aes(x = Range1,
                         y = jaccard.dist,
                         fill = Range1)) +
  geom_violin(color = 'black') +
  geom_boxplot(width = 0.1) +
  scale_fill_manual(values=c("lightblue", "grey")) +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        text = element_text(size=20),
        legend.position = "none",
        axis.ticks.y = element_blank())

jacc.plot <- ggpubr::ggarrange(btwn.plot, wtn.plot)
# ggsave('Fig2B.tiff',
#        plot = jacc.plot, device = 'tiff', scale = 2,
#        width = 12, height = 10, unit = 'cm')



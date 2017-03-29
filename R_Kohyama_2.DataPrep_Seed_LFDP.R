#' ---
#' title: "R_Kohyama_2.DataPrep_Seed_LFDP.R"
#' author: "Maria Wang"
#' date: "2017-03-27"
#' description: "This script preps seed rain data to be used in downstream analyses in other scripts"
#' ---

# ~ Loading required R packages ----
library(dplyr)
library(tidyr)

# ~ Set working directory ----
setwd("C:/Users/wang/Dropbox/LFDP/data analysis for Seth/Kohyama Hypothesis") #change this as needed
getwd()

# ~ Create folder for input data ----
dir.create("C:/Users/wang/Dropbox/LFDP/data analysis for Seth/Kohyama Hypothesis/input data/lfdp", recursive=TRUE)

# ~ Create folders to store output ----
dir.create("C:/Users/wang/Dropbox/LFDP/data analysis for Seth/Kohyama Hypothesis/figures.seed/")
dir.create("C:/Users/wang/Dropbox/LFDP/data analysis for Seth/Kohyama Hypothesis/intermediate data/errors/", recursive=TRUE)

# ~ Load required input files----
phenall <- read.csv("input data/lfdp/LFDP2_ElVerdePhenology_2006-2015.csv") #phenology basket data
vines.palms <- read.csv("input data/lfdp/lfdp_vines.palms.csv") #this is a list of species that are vines or palms or should not be included (e.g. UNKSPP, COFARA, Clusia sp, )
mult <- read.csv("input data/lfdp/lfdp_phenologyfruittoseedmultipliers.csv") # Fruit to Seed Multiplier
height <- read.csv("input data/lfdp/lfdp_tree_heights.csv")  # Max tree heights (m) by sp
smass <- read.csv("input data/lfdp/lfdp_seed_mass.csv") # Seed mass (g) by sp


#### Data Prep---- ####
str(phenall)
range(phenall$BASKET)
str(mult)

# ~ Rename species column to 'sp' to match other dataframes----
names(phenall) <- sub('SPECIES', 'sp', names(phenall)) 
names(mult) <- sub('SPP.CODE', 'sp', names(mult)) 

# ~ Omit vines and palms----
phenall <- phenall %>% anti_join(vines.palms) 

# ~ Convert date for easier handling----
phenall <- separate(phenall, DATE, into = c("month", "day", "year"), sep="/")
str(phenall)
phenall$month <- as.integer(phenall$month)
phenall$day <- as.integer(phenall$day)
phenall$year <- as.integer(phenall$year)
phenall$sp <- as.factor(phenall$sp)


#### ~ Separate into year-long intervals from April through March ####
phenall$interval.endyear <- 0
phenall[ (phenall$year == 2007 & phenall$month > 3) | (phenall$year == 2008 & phenall$month < 4), ]$interval.endyear <- 2008
phenall[ (phenall$year == 2008 & phenall$month > 3) | (phenall$year == 2009 & phenall$month < 4), ]$interval.endyear <- 2009
phenall[ (phenall$year == 2009 & phenall$month > 3) | (phenall$year == 2010 & phenall$month < 4), ]$interval.endyear <- 2010
phenall[ (phenall$year == 2010 & phenall$month > 3) | (phenall$year == 2011 & phenall$month < 4), ]$interval.endyear <- 2011
phenall[ (phenall$year == 2011 & phenall$month > 3) | (phenall$year == 2012 & phenall$month < 4), ]$interval.endyear <- 2012
phenall[ (phenall$year == 2012 & phenall$month > 3) | (phenall$year == 2013 & phenall$month < 4), ]$interval.endyear <- 2013
phenall[ (phenall$year == 2013 & phenall$month > 3) | (phenall$year == 2014 & phenall$month < 4), ]$interval.endyear <- 2014
phenall[ (phenall$year == 2014 & phenall$month > 3) | (phenall$year == 2015 & phenall$month < 4), ]$interval.endyear <- 2015

phenall <- phenall %>% filter(interval.endyear > 0) 
str(phenall)

#### ~ Select only seed/mature fruit data ####
# Phenology codes defined as:
# 0 = empty basket
# 1 = flower
# 2 = aborted fruit 
# 3 = immature fruit
# 4 = mature fruit --useful
# 5 = seed --useful
# 6 = pedicel of fruit 

phenall <- phenall %>% filter(CODE == 4 | CODE == 5) %>% #filter only mature fruit and seeds
  filter(is.na(NUMBER) == FALSE) #also exclude where the seeds/fruits were not counted
table(phenall$CODE)


#### ~ Number of fruits and seeds by species and year-interval ####
nrepro <- phenall %>% group_by(sp, CODE, interval.endyear) %>% summarize(sum=sum(NUMBER))


### ~ Join phenology and multiplier dataframes by species ####
phen <- phenall %>% left_join(mult, by="sp") %>% #join
  filter(is.na(MULTIPLIER)==FALSE) %>%          #omit species without multiplier
  rename(cohort = interval.endyear, census = CENSUS, code = CODE, plotno = BASKET ) %>%
  arrange(sp, cohort, census, code, plotno) #sort data by species and date
head(phen)

n_distinct(phen$sp) #N=49
unique((phen$sp))


#### ~ Filter out species missing multipliers ####
spp.no.mult <- phenall %>% left_join(mult, by="sp") %>% #join
  filter(is.na(MULTIPLIER)==TRUE) %>% distinct(sp)

spp.no.mult2 <- inner_join(spp.no.mult, height) %>% distinct(sp) #filter to only species with tree heights

write.csv(spp.no.mult2, "intermediate data/errors/species.without.multiplier.csv", row.names=FALSE)


#### ~ Calculate number of seeds by multiplying no. of fruits (code =4) with multiplier ####
phen$nseeds <- 0
phen$nseeds <- ifelse(phen$code==4, phen$NUMBER*phen$MULTIPLIER, phen$NUMBER)
head(phen)


#### ~ Filter to only include species with 10 or more seeds per year ####
spplist2 <- phen %>% group_by(sp, cohort) %>% summarize(tot.seeds.year = sum(nseeds)) %>% filter(tot.seeds.year > 9) 
seeds.n10 <- left_join(spplist2, phen)
n_distinct(seeds.n10$sp) #DOWN TO 38 spp if include more than 9
unique((seeds.n10$sp))

write.csv(seeds.n10, "intermediate data/lfdp_seedsprepped.n10.csv", row.names=FALSE)

#
#
#

######## ---- SEED TO SEEDLING (STS) ESTABLISHMENT ---- #######
# following Muscarella et al. 2013,
# for the GLMM, we need:
# response variable = no. of seedlings recruited in each seedling plot
# offset = no. of seeds recorded in each associated seed basket for corresponding seedling cohort year
# fixed effects = seed mass, spp height
# random effects = sp, cohort, plotID

# We have 3 x 1m2 sling plots per 1 x 0.5-m2 seed trap

# we calculated mean per-seed success of the seed-to-seedling (STS) transition for each species, henceforth referred to as "STS transition", as the number of seedling recruits found in each of the three 1 m2 seedling plots divided by six times the number of conspecific seeds observed in the associated 0.5 m2 basket for each corresponding year. 


#### ~ Load seedling and seed data ####
sling.survival.y1 <- read.csv("intermediate data/lfdp_seedling_forSTS.csv") #generated by R_Kohyama_1.DataPrep_Seedling_LFDP.R

seeds.n10 <- read.csv("intermediate data/lfdp_seedsprepped.n10.csv") #generated by code above

str(sling.survival.y1)
str(seeds.n10)

# ~ Group seed and seedling data by sp, cohort, plotno
sling.survival.y1b <- sling.survival.y1 %>% group_by(sp, cohort, plotno) %>% summarize(n.new.sling=n())
seeds.n10b <- seeds.n10 %>% group_by(sp, cohort, plotno) %>% summarize(nseeds=sum(nseeds))

# ~ Join seed and seedling dataframes
ss <- left_join(seeds.n10b, sling.survival.y1b) #we are calculating seed survival, so the first dataframe has all that we want to join
head(ss)

table(ss$cohort, ss$plotno) #check to see if each plot had data for each cohort. Yes they do! 

# ~ when no. of new slings is NA, convert to zero
ss[is.na(ss$n.new.sling), ]$n.new.sling <- 0


#### ~ Join STS, tree height and seed mass dataframes ####
sshm <- left_join(ss, height) %>% left_join(smass)
names(sshm)
sshm <- sshm %>% dplyr::select(sp, max.ht, seed.mass = seed.mass_in_g, cohort, plotno, nseeds, n.new.sling)
n_distinct(sshm$sp) #n sp = 38

#### ~ Check which species are missing max.ht and seed mass ####
sp.no.height <- sshm %>% filter( is.na(max.ht) == TRUE)  %>% distinct(sp)
sp.no.smass <- sshm %>% filter( is.na(seed.mass) == TRUE)  %>% distinct(sp)
write.csv(sp.no.smass, "intermediate data/errors/species.need.seedmass-seed.csv", row.names=FALSE)

#### ~ Filter out NAs ####
sshm <- sshm %>% filter(!is.na(max.ht) == TRUE)
sshm <- sshm %>% filter(!is.na(seed.mass) == TRUE)

#### ~ Convert categories to factor ####
sshm$sp <- as.factor(sshm$sp)
sshm$cohort <- as.factor(sshm$cohort)
sshm$plotno <- as.factor(sshm$plotno)


#### ~ Create tree height categories ####
# sshm$htclass <- factor(NA, c("understory", "canopy"))
# sshm[sshm$max.ht <= 15, ]$htclass <- "understory"
# sshm[sshm$max.ht > 15, ]$htclass <- "canopy"
# str(sshm$htclass)

#### ~ Standardize/Center max.ht and seed.mass ####
sshm$max.ht.z <- (sshm$max.ht - mean(sshm$max.ht) )/(2*sd(sshm$max.ht))
sshm$seed.mass.z <- (sshm$seed.mass - mean(sshm$seed.mass, na.rm=TRUE) )/(2*sd(sshm$seed.mass, na.rm=TRUE))

#### ~ Create data for GLMM ####
check <- sshm[is.na(sshm$n.new.sling)==TRUE | is.na(sshm$nseeds)==TRUE, ] #checking for missing data
sshm2 <- sshm[is.na(sshm$n.new.sling)==FALSE, ]  #omit NAs

#### ~ When n.new.sling > nseeds, set total number of seeds divided by 0.5m2 * 3 1-m2 sling plots = total number of new seedling recruits (following Muscarella et al. 2013)
sshm2$nseeds.x6 <- sshm2$nseeds/0.5*3 #need to multiply # seeds by 6 to make it equal area with seedling plots
sshm2$nseeds.x6 <- ifelse(sshm2$n.new.sling > sshm2$nseeds.x6, sshm2$n.new.sling, sshm2$nseeds.x6) #if n.new.sling > nseeds.x6, replace value of nseeds.x6 with value of n.new.sling so they are equal.

sshm2$failures <- sshm2$nseeds.x6 - sshm2$n.new.sling
sshm2$sts <- sshm2$n.new.sling/sshm2$nseeds.x6 #sts = mean per-seed-success
head(sshm2)

write.csv(sshm2, "intermediate data/ldfp_sts_forglm.csv", row.names=FALSE)

#############(not needed for GLM)################
#### Mean per-seed-success by species #### 
stsbyspp <- sshm2 %>% group_by(sp) %>% summarize(n.new.sling = sum(n.new.sling), nseeds.x6 = sum(nseeds.x6), sts=sum(n.new.sling)/sum(nseeds.x6))


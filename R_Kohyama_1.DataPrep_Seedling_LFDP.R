#' ---
#' title: "R_Kohyama_1.DataPrep_Seedling_LFDP.R"
#' author: "Maria Wang"
#' date: "2017-03-27"
#' description: "This script preps seedling data to be used in downstream analyses in other scripts"
#' ---

# ~ Load required R packages ----
library(dplyr)
library(tidyr)

# ~ Set working directory ----
setwd("C:/Users/wang/Dropbox/LFDP/data analysis for Seth/Kohyama Hypothesis") #change this as needed
getwd()

# ~ Create folder for input data ----
dir.create("C:/Users/wang/Dropbox/LFDP/data analysis for Seth/Kohyama Hypothesis/input data/lfdp", recursive=TRUE)

# ~ Create folders to store output ----
dir.create("C:/Users/wang/Dropbox/LFDP/data analysis for Seth/Kohyama Hypothesis/figures.seedling/")
dir.create("C:/Users/wang/Dropbox/LFDP/data analysis for Seth/Kohyama Hypothesis/intermediate data/errors/", recursive=TRUE)


# ~ Load required input files----
sling0 <-read.csv("input data/lfdp/LFDP_PhenologySeedlings_2007-2015MW.csv") # LFDP seedling census data
vines.palms <- read.csv("input data/lfdp/lfdp_vines.palms.csv") # list of species that are vines or palms or should not be included (e.g. UNKSPP, COFARA, Clusia sp, )

str(sling0)

#### Data Prepping ####
# Create a new column with unique ID number for each seedling
sling0$uniqueid <- 1:nrow(sling0)

# Exclude seedlings from 2007 Census (because unknown age)
sling0 <- filter(sling0, is.na(CENSUS1_2007) == TRUE) 

# Subset data to include only identifiers and survival status 
colnames(sling0)
head(sling0)
sling <- select(sling0, uniqueid, plotno = PLOTNUMBER, subplot = SUBPLOT, plotid = PLOTID, tag = TAG, sp = FINALSPP2007_2014, SUMMARYSTATUS2014, A.D.N.NF.Q.R2008, A.D.N.NF.Q.R2009, A.D.N.NF.Q.R2010, A.D.N.NF.Q.R2011, A.D.N.NF.Q.R2012, A.D.N.NF.Q.R2013, A.D.N.NF.Q.R2014, A.D.NF.Q.R2015, STATUS2008, STATUS2009, STATUS2010, STATUS2011, STATUS2012, STATUS2013, STATUS2014, STATUS2015)

table(sling$SUMMARYSTATUS2014)

# Exclude R and Q
sling <- filter(sling, !(SUMMARYSTATUS2014 == "R" | SUMMARYSTATUS2014 == "Q" ) )
sling <- filter(sling, !(STATUS2008 == "R" | STATUS2008 == "Q" | STATUS2009 == "R" | STATUS2009 == "Q" | STATUS2010 == "R" | STATUS2010 == "Q" | STATUS2011 == "R" | STATUS2011 == "Q" | STATUS2012 == "R" | STATUS2012 == "Q" | STATUS2013 == "R" | STATUS2013 == "Q" | STATUS2014 == "R" | STATUS2014 == "Q" | STATUS2015 == "R" | STATUS2015 == "Q") )
sling <- droplevels(sling)

# omit vines and palms 
sling <- sling %>% anti_join(vines.palms)
n_distinct(sling$sp) #N=57
unique((sling$sp))

#### ~ No. of New Seedlings by Spp and Cohort (FOR SEEDLING SURVIVAL!) ####
sling$cohort <- 0
# assign cohort to year X if NEW ('N') AND ALIVE ('A') in that year
sling[(sling$A.D.N.NF.Q.R2008 == 'N' & sling$STATUS2008=='A'), ]$cohort <- 2008
sling[(sling$A.D.N.NF.Q.R2009 == 'N' & sling$STATUS2009=='A'), ]$cohort <- 2009
sling[(sling$A.D.N.NF.Q.R2010 == 'N' & sling$STATUS2010=='A'), ]$cohort <- 2010
sling[(sling$A.D.N.NF.Q.R2011 == 'N' & sling$STATUS2011=='A'), ]$cohort <- 2011
sling[(sling$A.D.N.NF.Q.R2012 == 'N' & sling$STATUS2012=='A'), ]$cohort <- 2012
sling[(sling$A.D.N.NF.Q.R2013 == 'N' & sling$STATUS2013=='A'), ]$cohort <- 2013
sling[(sling$A.D.N.NF.Q.R2014 == 'N' & sling$STATUS2014=='A'), ]$cohort <- 2014
sling[(sling$A.D.NF.Q.R2015 == 'N' & sling$STATUS2015=='A'), ]$cohort <- 2015

sling <- sling %>% filter(cohort > 0) #omit seedlings unassigned to cohorts

#### ~ 1st Year Seedling Survival (% y-1) for model #######
#King 2006: First-year seedling survival (% y-1) was estimated as the proportion
#of new recruits that survived until the next annual census for species with 10 or more seedling recruits.

# For each recruit, I assign a '1' to the new recruits that survived until the next annual census = creating binary data for glmm
sling$live.in.y1 <- 0
sling[(sling$cohort == 2008 & sling$STATUS2009 == 'A'), ]$live.in.y1 <- 1
sling[(sling$cohort == 2009 & sling$STATUS2010 == 'A'), ]$live.in.y1 <- 1
sling[(sling$cohort == 2010 & sling$STATUS2011 == 'A'), ]$live.in.y1 <- 1
sling[(sling$cohort == 2011 & sling$STATUS2012 == 'A'), ]$live.in.y1 <- 1
sling[(sling$cohort == 2012 & sling$STATUS2013 == 'A'), ]$live.in.y1 <- 1
sling[(sling$cohort == 2013 & sling$STATUS2014 == 'A'), ]$live.in.y1 <- 1
sling[(sling$cohort == 2014 & sling$STATUS2015 == 'A'), ]$live.in.y1 <- 1

str(sling)

# ~ Filter to only include species with 10 or more new seedling recruits per year ####
spplist2 <- sling %>% group_by(sp, cohort) %>% summarize(tot.new.sling = n()) %>% filter(tot.new.sling > 9) 
sling.n10 <- left_join(spplist2, sling)
n_distinct(sling.n10$sp) #DOWN TO 23 spp
unique((sling.n10$sp))


#### ~ Save tables as csv for use in further analyses ####

# This csv will be used in R_Kohyama_2.DataPrep_Seed_LFDP.R
write.csv(sling, "intermediate data/lfdp_seedling_forSTS.csv", row.names=FALSE) 

# This csv will be used in R_Kohyama_4.ModelTesting_Seedling_LFDP.R
write.csv(sling.n10, "intermediate data/lfdp_seedling_excldlessthan10recruits.csv", row.names=FALSE) 


##### Format Data for GLM Model Testing ####
sling.n10 <- read.csv("intermediate data/lfdp_seedling_excldlessthan10recruits.csv")
names(sling.n10)
str(height)
str(smass)

#### ~ Join seedling survival, tree height and seed mass data ####
height <- read.csv("input data/lfdp/lfdp_tree_heights.csv") # Max tree heights (m) by sp
smass <- read.csv("input data/lfdp/lfdp_seed_mass.csv") # Seed mass (g) by sp

slinghm <- left_join(sling.n10, height) %>% left_join(smass)
names(slinghm)
slinghm <- slinghm %>% dplyr::select(uniqueid, sp, max.ht, seed.mass = seed.mass_in_g, cohort, plotid, plotno, subplot, live.in.y1)
n_distinct(slinghm$sp) #n sp = 23

#### ~ Check which species are missing max.ht and seed mass ####
sp.no.height <- slinghm %>% filter( is.na(max.ht) == TRUE)  %>% distinct(sp)
sp.no.smass <- slinghm %>% filter( is.na(seed.mass) == TRUE)  %>% distinct(sp)
write.csv(sp.no.smass, "intermediate data/errors/species.need.seedmass-sling.csv", row.names=FALSE)

#### ~ Filter out NAs ####
slinghm <- slinghm %>% filter(!is.na(max.ht) == TRUE)
slinghm <- slinghm %>% filter(!is.na(seed.mass) == TRUE) #3 sp
n_distinct(slinghm$sp) #n sp = down to 20

#### ~ Convert categories to factor ####
slinghm$sp <- as.factor(slinghm$sp)
slinghm$cohort <- as.factor(slinghm$cohort)
slinghm$plotid <- as.factor(slinghm$plotid)
slinghm$plotno <- as.factor(slinghm$plotno)

#### ~ Create tree height categories ####
# slinghm$htclass <- factor(NA, c("understory", "canopy"))
# slinghm[slinghm$max.ht <= 15, ]$htclass <- "understory"
# slinghm[slinghm$max.ht > 15, ]$htclass <- "canopy"
# str(slinghm$htclass)

#### ~ Standardize/Center max.ht and seed.mass ####
slinghm$max.ht.z <- (slinghm$max.ht - mean(slinghm$max.ht) )/sd(slinghm$max.ht)
slinghm$seed.mass.z <- (slinghm$seed.mass - mean(slinghm$seed.mass, na.rm=TRUE) )/sd(slinghm$seed.mass, na.rm=TRUE)

check <- slinghm[is.na(slinghm$live.in.y1)==TRUE, ] #checking to see if there is missing data
slinghm <- slinghm[is.na(slinghm$live.in.y1)==FALSE, ]  #omit NAs

#### ~ Create CSV for seedling survival GLM ####
write.csv(slinghm, "intermediate data/lfdp_seedlinghm_forglm.csv", row.names=FALSE)

######### Exploratory Stuff (not required for analyses) ----- #########
#### ~ Mean Seedling Survival by species ####
slinghm <- read.csv("intermediate data/lfdp_seedlinghm_forglm.csv")
names(slinghm)
slingsurvival.y1.sp <- slinghm %>% group_by(sp, max.ht, seed.mass, max.ht.z, seed.mass.z) %>% 
  summarize(n.new.sling=n(), n.live.y1=sum(live.in.y1)) %>%
  mutate(slingsurvival.y1 = n.live.y1/n.new.sling)

#### ~ Mean Seedling Survival by species, cohort, plot ####
slingsurvival.y1.all <- slinghm %>% group_by(sp, max.ht, seed.mass, max.ht.z, seed.mass.z, cohort, plotid, plotno, subplot) %>% 
  summarize(n.new.sling=n(), n.live.y1=sum(live.in.y1)) %>%
  mutate(slingsurvival.y1 = n.live.y1/n.new.sling)

#write.csv(slingsurvival.y1.all, "intermediate data/lfdp_seedlinghm_byrandeffects.csv", row.names=FALSE)

#### ~ Exploratory plots ####
plot(slingsurvival.y1 ~ max.ht, data=slingsurvival.y1.sp)
plot(slingsurvival.y1 ~ seed.mass, data=slingsurvival.y1.sp)
plot(log1p(slingsurvival.y1) ~ max.ht, data=slingsurvival.y1.sp)
plot(log1p(slingsurvival.y1) ~ log1p(seed.mass), data=slingsurvival.y1.sp)
plot(log10(slingsurvival.y1) ~ log10(seed.mass), data=slingsurvival.y1.sp) 
plot(slingsurvival.y1 ~ seed.mass, data=slingsurvival.y1.sp, log="xy")
hist(slingsurvival.y1.sp$slingsurvival.y1, breaks=50)

#### ~ Exploratory plots for averaged data ####
slingsurvival.y1.all <- read.csv("intermediate data/slinghm_byrandeffects.csv")
plot(slingsurvival.y1 ~ max.ht, data=slingsurvival.y1.all)
plot(slingsurvival.y1 ~ seed.mass, data=slingsurvival.y1.all)
plot(log1p(slingsurvival.y1) ~ max.ht, data=slingsurvival.y1.all)
plot(log1p(slingsurvival.y1) ~ log1p(seed.mass), data=slingsurvival.y1.all)
plot(log10(slingsurvival.y1) ~ log10(seed.mass), data=slingsurvival.y1.all) 
plot(slingsurvival.y1 ~ seed.mass, data=slingsurvival.y1.all, log="xy")
hist(slingsurvival.y1.all$slingsurvival.y1, breaks=50)


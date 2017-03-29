#### R_Kohyama_3.DataPrep_SaplingRecruitment_LFDP.R ####
#### Kohyama Hypothesis Project
#### Code by Maria Wang, Summer/Fall 2016

#### Load packages ####
library(dplyr)
library(tidyr)
library(date)
attach('C:/Users/wang/Dropbox/LFDP/data analysis for Seth/CTFSRPackage.rdata')
ls(2)

#### Set working directory ####
setwd("C:/Users/wang/Dropbox/LFDP/data analysis for Seth/Kohyama Hypothesis/")
getwd()

#### Load tree census data ####
load('input data/lfdp/lfdp_R.tables/luquillo.stem1.Rdata') 
load('input data/lfdp/lfdp_R.tables/luquillo.stem2.Rdata') 
load('input data/lfdp/lfdp_R.tables/luquillo.stem3.Rdata') 
load('input data/lfdp/lfdp_R.tables/luquillo.stem4.Rdata') 
load('input data/lfdp/lfdp_R.tables/luquillo.stem5.Rdata') 
str(luquillo.stem1)
unique(luquillo.stem1$sp)

#### Calculate recruitment using CTFS R package function ####
rec1 <- recruitment.eachspp(luquillo.stem1, luquillo.stem2, mindbh=15) %>% as.data.frame 
rec2 <- recruitment.eachspp(luquillo.stem2, luquillo.stem3, mindbh=15) %>% as.data.frame 
rec3 <- recruitment.eachspp(luquillo.stem3, luquillo.stem4, mindbh=15) %>% as.data.frame 
rec4 <- recruitment.eachspp(luquillo.stem4, luquillo.stem5, mindbh=15) %>% as.data.frame 

## Rename column names
rec1$sp <- row.names(rec1)
rec2$sp <- row.names(rec2)
rec3$sp <- row.names(rec3)
rec4$sp <- row.names(rec4)

### Combine dataframes ####
rec.all <- bind_rows(rec1, rec2, rec3, rec4)
rec.all$interval <- rep(1:4, each=nrow(rec1)) #create column indicating interval no.
names(rec.all)

## Clean up 
rec.all2 <- rec.all[is.na(rec.all$rate) == FALSE, ] #omit 121 NAs
range(rec.all2$rate) #Inf will cause issues when log-transforming later
rec.all2 <- rec.all2[rec.all2$rate < Inf, ] #omit Inf

## Only incld spp where number of individuals alive and >15mm dbh in census 2 is at least 10 
rec.all2 <- rec.all2[rec.all2$N2 > 9, ] 
rec.all0 <- rec.all2[rec.all2$rate == 0, ] #20 zeros

write.csv(rec.all2, "output.sapling/lfdp_recruit_spp.csv", row.names=FALSE)

#### Preparing data for model testing ####
rec.all2 <- read.csv("output.sapling/lfdp_recruit_spp.csv")
treeht <- read.csv("input data/lfdp/lfdp_tree_heights.csv") 
str(rec.all2)
str(treeht)
unique(rec.all2$sp)
unique(treeht$sp)

#omit vines and palms
vines.palms <- read.csv("input data/lfdp/vines.palms.csv") #this is a list of species that are vines or palms or should not be included (e.g. UNKSPP, COFARA, Clusia sp, )
vines.palms <- read.csv("https://dl.dropboxusercontent.com/u/43806775/input%20data/vines.palms.csv")
rec.all2 <- rec.all2 %>% anti_join(vines.palms)
n_distinct(rec.all2$sp) #N=87
unique((rec.all2$sp))

#check what tree heights are missing
need_heights <- anti_join(rec.all2, treeht) %>% distinct(sp)

write.csv(need_heights, "output.sapling/need_heights_lfdp.csv")

# join dataframes
recruit <- inner_join(rec.all2, treeht) #only keep sp with max heights
str(recruit)
head(recruit)

# write csv ####
write.csv(recruit, "output.sapling/lfdp_recruit_lme.csv", row.names=FALSE)

recruit[is.na(recruit$rate) == TRUE, ] #checking
range(recruit$rate)
